function Top_mfse_Multi_field_2D(nptx,npty,refine,volfra,corlencoe,nums_x,nums_y)
%% MATERIAL AND FIELD PROPERTIES
E0 = 200.0e5; Emin = 1e-9*E0; nu = 0.3; ptdist = 1;
elelen = ptdist/refine; nelx = refine*nptx; nely = refine*npty;
tolne = nelx*nely; tolnd = (nelx+1)*(nely+1); tolnpt = nptx*npty;
tolvol = tolne*elelen^2;
fprintf([' Number of material-field  points:%10i \n'...
    ' Number of finite elements:%10i\n'],tolnpt,tolne);

%% PREPARE FINITE ELEMENT ANALYSIS
nodenrs = reshape(1:tolnd,1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,tolne,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1], tolne,1);
iK = reshape(kron(edofMat,ones(8,1))', 64*tolne,1);
jK = reshape(kron(edofMat,ones(1,8))', 64*tolne,1);
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE0 = 1/(1-nu^2)/24*([A11 A12; A12' A11]+nu*[B11 B12; B12' B11]);
% SETUP BOUNDARY CONDITIONS (HALF MBB-BEAM)
F = sparse(2, 1, -1000, 2*tolnd, 1); U = zeros(2*tolnd,1);
fixeddofs = [1:2:2*(nely + 1),2*(nelx + 1)*(nely + 1),2*(nelx)*(nely + 1)];
freedofs = setdiff(1:2*tolnd,fixeddofs);

%% PREPARE MATERIAL FIELD SERIES EXPANSION
[eIntopMat,ptIntopMat] = MFSE2D(nptx/nums_x,npty/nums_y,refine,corlencoe);

%% INITIALIZE DESIGN VARIABLES
beta = 0.5; penal = 3;
x = (-log(1/volfra-1)/beta)*ones(1,tolnpt/(nums_x*nums_y))/ptIntopMat;
x = x'; neig1 = length(x); 
x = repmat(x,nums_x*nums_y,1); neig2 = length(x); n = neig2; m = 1;

%% INITIALIZE MMA OPTIMIZER
loop = 0; obj = 0.;
change = 1.; ichange = 1;
xmin = -1000*ones(n,1); xmax = 1000*ones(n,1);
low = xmin; upp = xmax;
xold1 = x;  xold2 = x; clf;
cc = 10000*ones(m,1); d = zeros(m,1); a0 = 1; a = zeros(m,1);
Obj = []; Volf = [];
[Xe, Ye] = meshgrid((0.5:1:nelx)*elelen, (nely-0.5:-1:0.5)*elelen);

%% START ITERATION
while (change >= 0.005 || beta < 20)
    loop = loop + 1; objold = obj;
    %% MATERIAL FIELDS IN THE SUB-DOMAIN
    x_design = cell(nums_y,nums_x); ePhi = cell(nums_y,nums_x); 
    ePhiProj = cell(nums_y,nums_x); edproj = cell(nums_y,nums_x);
    temp = 0;
    for i = 1 : nums_y
        for j = 1 : nums_x
            x_design{i,j} = x(temp + 1:temp + neig1);
            ePhi{i,j} = eIntopMat'*x_design{i,j};
            [ePhiProj{i,j},edproj{i,j}] = threshold(ePhi{i,j},beta);
            temp = temp + neig1;
        end
    end
    
    %% ASSEMBLE THE MATERIAL FIELDS IN THE DESIGN DOMAIN
    temp_y = 0; phi = zeros(nely,nelx); plot_ePhi = zeros(nely,nelx);
    for i = 1 : nums_y
        temp_x = 0;
        for j = 1 : nums_x
            phi_temp = reshape(ePhiProj{i,j},nely/nums_y,nelx/nums_x);
            phi(temp_y + 1 : temp_y + nely/nums_y, temp_x+1:temp_x + nelx/nums_x) = phi_temp;
            plot_ePhi(temp_y + 1 : temp_y + nely/nums_y, temp_x+1:temp_x + nelx/nums_x) = reshape(ePhi{i,j},nely/nums_y,nelx/nums_x);
            temp_x = temp_x + nelx/nums_x;
        end
        temp_y = temp_y + nely/nums_y;
    end
    
    %% FE-ANALYSIS
    sK = reshape(KE0(:)*(Emin + (E0 - Emin)*phi(:)'.^penal), 64*tolne, 1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    obj = F'*U;
    ce_temp = reshape(sum((U(edofMat)*KE0).*U(edofMat), 2),nely,nelx);
    ce = cell(nums_y,nums_x); dcdx = cell(nums_y,nums_x); voldgdx = cell(nums_y,nums_x);
    temp_y = 0; dc = []; dv = []; x = [];
    for i = 1 : nums_y
        temp_x = 0;
        for j = 1 : nums_x
            ce{i,j} = ce_temp(temp_y + 1 : temp_y + nely/nums_y, temp_x+1:temp_x + nelx/nums_x);
            ePhiProj_temp = ePhiProj{i,j};
            dcdx{i,j} = eIntopMat*(-penal*(E0 - Emin)*(ePhiProj_temp(:)).^(penal-1).*ce{i,j}(:).*edproj{i,j}(:));
            voldgdx{i,j} = eIntopMat*(edproj{i,j}(:)*elelen^2);
            temp_x = temp_x + nelx/nums_x;
            dc = cat(1,dc,dcdx{i,j}); dv = cat(1,dv,voldgdx{i,j}); x = cat(1,x,x_design{i,j});
        end
        temp_y = temp_y + nely/nums_y;  
    end
    vol = sum((phi(:))*elelen^2);
    
    %% UPDATE BY THE MMA OPTIMIZER
    fval = zeros(m, 1); fval(1) = 100*(vol/tolvol-volfra);
    dfdx = zeros(m, n); dfdx(1,:) = 100*dv/tolvol;   
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,loop,x,xmin,xmax,xold1,xold2,obj,dc,fval,dfdx,low,upp,a0,a,cc,d);
    xold2 = xold1; xold1 = x; x = xmma;
    
    %% TUNE PROJECTION PARAMETER
    change = abs(obj-objold)/obj;
    if change < 0.005 && loop > 30
        ichange = ichange+1;
    else
        ichange = 1;
    end
    if mod(ichange,3) == 0
        beta = min(beta * 1.1,20);
    end
    
    %% PRINT RESULTS AND PLOT DENSITIES
    fprintf([' It.:%5i Obj.:%9.4f Vol:%7.4f numdesvars :%5i' ...
        ' beta:%5.1f ch.:%6.3f\n'],...
        loop,obj,vol/tolvol,neig2,beta,change);  
    figure(1); clf;
    displayx = zeros(nely, 2*nelx);
    displayx(:, 1:nelx) = flip(reshape(phi, nely, nelx),2);
    displayx(:, nelx+1:end) = displayx(:, nelx:-1:1);
    colormap(gray); clims=[-1 0]; imagesc(-displayx,clims); 
    axis equal; axis tight; title('Elemental density distribution');
    set(gca,'XTick',[0 1e5]);set(gca,'YTick',[0 1e5]);
    figure(2); clf;
    plot_ePhi = reshape(plot_ePhi,size(Xe));
    contourf([Xe nptx*ptdist+Xe],[Ye Ye],[fliplr(plot_ePhi) plot_ePhi],[0 0]);
    colormap([0 0 0; 0 0 1;1 0 0; 0 1 0; 1 1 1]);
    title('Material-field contour of the optimized design');
    axis equal; axis tight; set(gca,'XTick',[]); set(gca,'YTick',[]);
    figure(3); clf;
    Obj = cat(2,Obj,obj); Volf = cat(2,Volf,vol/tolvol);
    plotConvergence(Obj,Volf);
end
end

function [eIntopMat,ptIntopMat] = MFSE2D(nptx,npty,refine,corlencoe)
ptdist = 1; corlen = corlencoe*min(nptx,npty)*ptdist;
elelen = ptdist/refine; nelx = refine*nptx; nely = refine*npty;
tolne = nelx*nely;  tolnpt = nptx*npty;
%% BUILD CORRELATION MATRIX
[Xpt, Ypt] = meshgrid((0.5:1:nptx)*ptdist, (npty-0.5:-1:0.5)*ptdist);
Xpt = Xpt(:); Ypt = Ypt(:);
corMat = zeros(tolnpt,tolnpt);
for i = 1:size(corMat,1)
    for j = i+1:size(corMat,2)
        corMat(i,j) = exp(-(((Xpt(j)-Xpt(i))^2+(Ypt(j)-Ypt(i))^2)/corlen^2));
    end
end
corMat = corMat+corMat';
for i = 1:size(corMat, 1)
    corMat(i,i) = 1;
end
%% DO SERIES EXPANSION OF THE MATERIAL FIELD
if size(corMat,1) < 1e4
    [eigfunMat, eigvalMat] = eig(corMat);
else
    [eigfunMat, eigvalMat] = eigs(corMat,1500);
end
eigvalVec = diag(eigvalMat);
[eigvalVec, eigsortind] = sort(eigvalVec, 'descend');
neig = 0; tmpsum = 0.;
while tmpsum < (1-1e-4)*sum(abs(eigvalVec))
    neig = neig + 1;
    tmpsum = tmpsum + eigvalVec(neig);
end
EXPANMat = sparse(1:neig, 1:neig, eigvalVec(1:neig).^(-1/2), neig, neig)...
    *eigfunMat(:,eigsortind(1:neig))'; clear eigfunMat;
%% COMPUTE PHI ON ELEMENTS AND MATERIAL-FIELD POINTS
[Xe, Ye] = meshgrid((0.5:1:nelx)*elelen, (nely-0.5:-1:0.5)*elelen);
Xe = Xe(:); Ye = Ye(:);
eIntopMat = zeros(neig, tolne);
grsize = min(round(tolnpt/20), tolne); ngr = ceil(tolne/grsize);
for igr = 1:ngr
    eind = (igr-1)*grsize+1:min(igr*grsize, tolne);
    Xe_sub = Xe(eind); Ye_sub = Ye(eind);
    eptvals = exp(-(((repmat(Xpt',length(eind),1)-repmat(Xe_sub, 1, tolnpt)).^2 ...
        +(repmat(Ypt',length(eind),1)-repmat(Ye_sub, 1, tolnpt)).^2)/corlen^2))';
    eptvals(abs(eptvals)<1e-9) = 0;
    eIntopMat(:,eind) = EXPANMat*eptvals;
end
ptIntopMat = EXPANMat*corMat'; clear corMat;
end

function [ePhiProj, edproj] = threshold(ePhi, beta)
%% SIGMOID PROJECTION
ePhiProj = 1./(1+exp(-beta*ePhi));
edproj = beta*ePhiProj.*(1-ePhiProj);
end

function plotConvergence(values , volfrac)
iter = 1:length(values);
%% COMPLIANCE
axes1 = gca;
yyaxis(axes1,'left');
plot(iter, values,'b-','LineWidth',1.5);
ylabel('Structural compliance','FontSize',14,'FontName','Times New Roman','Color',[0 0 1]);
set(axes1,'YColor','b','FontSize',14,'FontName','Times New Roman');

%% VOLUME FRACTION
yyaxis(axes1,'right');
plot(iter,volfrac,'r-.','LineWidth',1.5);hold on;
set(axes1,'ylim',[0 1]);set(axes1,'ytick',0:.1:1);
ylabel('Volume constraint','FontSize',14,'FontName','Times New Roman','Color',[1 0 0]);
set(axes1,'YColor','r','FontSize',14,'FontName','Times New Roman');

%% AXIS
xlabel('Number of iterations','FontSize',14,'FontName','Times New Roman');
drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is supplementary to the corresponding paper:                    %
% A multi-regional MFSE method for large-scale structures                  %
% with arbitrary design domain                                             %
% Zhaoyou Sun, Tingxi Yuan, Wenbo Liu, Jiaqi He, Tiejun Sui, Yangjun Luo   %
% This code is based on the published educational paper                    %
% A Matlab Code for the Material-Field Series-Expansion Topology           %
% Optimization Method, by Liu et al., Front. Mech. Eng. (2021)             %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guaranty that the code is      %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
