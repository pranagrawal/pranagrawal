
function Bsurfexample

clear all
clc
close all

% B surface level set function, zero level set is hole

set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontSize',10)
colormap parula


%domain
xori = 0.;
yori = 0.;
xdom = 1.;
ydom = 1.;

% number of control points in each direction
nvarx = 15;
nvary = 15;

% number of patches
nexv = nvarx - 3;
neyv = nvary - 3;

% number of variables
nvar = nvarx*nvary;

%size of patch Bsurface elements
dxv=(xdom)/nexv;
dyv=(ydom)/neyv;

%Control points coordinates location
xplot=linspace(xori-dxv,xori+xdom+dxv,nvarx);
yplot=linspace(yori-dyv,yori+ydom+dyv,nvary);

%mesh of control points in x and y
[Xcs,Ycs]=meshgrid(xplot,yplot);

% center of the domain
xc = xori + .5*xdom;
yc = yori + .5*ydom;

% star like (example)
data.myFigname = 'FigStar';
theta = atan((Ycs-yc)./(Xcs-xc));
radius = 0.4 + 0.06*sin(6*theta);
dvarM=5*(radius - sqrt((Xcs-xc).^2 + (Ycs-yc).^2));

% distance between points to get zero contour
dl = 0.01;

%% pack
data.xori = xori;
data.yori = yori;
data.xdom = xdom;
data.ydom = ydom;
data.nexv = nexv;
data.neyv = neyv;
data.nvar = nvar;
data.dl = dl; % distance of separation between points
data.Xcs = Xcs;
data.Ycs = Ycs;
data.myFontSizeLabel = 12;
data.sizecmf = [7 7];

% Matrix to vector form of design variables
dvar= dvarM(:);

% ploting Bsurface
plotBsurf(data,dvarM);

% find zero contour
points = findedge(data, dvar)
%specifying the mesh parameters

lc = 0.01; %defien the leact count of the mesh size

% 2D mesh algorithm 
% 1: MeshAdapt, 
% 2: Automatic,
% 3: Initial mesh only, 
% 5: Delaunay,
% 6: Frontal-Delaunay, 
% 7: BAMG, 
% 8: Frontal-Delaunay for Quads, 
% 9: Packing of Parallelograms 
% Default value: 6
malgo = 8;
meshplot(points,lc,malgo);

end


function points = findedge(data, dvar)

% distance of separation between points
dl = data.dl;

% maximum number of iterations
nmax = 50;
tolmin = 1.d-8;

% centroid
xm = data.xori+.5*data.xdom;
ym = data.yori+.5*data.ydom;

% point a at center
xa = xm;
ya = ym;

% point b at right side
xb = data.xori+data.xdom;
yb = ym;


% starting point for edge traverse
[xo,dphip] = startingp(xa,ya,xb,yb, data, tolmin, dvar);


% lambda
lambda = 0.;
npoi = 1;
npoimax = 1000;
xinit = xo;

% store points
points = xo;
bcontinue = true;

while  npoi < npoimax && bcontinue
    % direction perpendicular
    mydire = [dphip(2) -dphip(1)];
    
    % guess point
    xp = xo + dl*mydire/norm(mydire);
    
    %find new point
    [xp,dphip,lambda]=findnewpoint(data,nmax,tolmin,xo,xp,lambda, dl, dvar);
    
    % distance to initial point
    
    
    %else
    % store
    points = [points;xp];
    % update
    xo = xp;
    npoi = npoi+1;
    
    if npoi>2 && ( norm(xp - xinit)<dl  )
        bcontinue = false;
    end
end

% plottin zero

figure()
plot(points(:,1),points(:,2),'ro','MarkerSize',1,'MarkerFaceColor','r')
hold on
%text(points(1,1),points(1,2),'${\bf Y}_1$','FontSize',data.myFontSizeLabel)
%annotation('\textarrow',points(1,1),points(1,2),'String','${\bf Y}_1$');
xlabel('$x_1$')
ylabel('$x_2$')
axis tight
axis equal
xlim([data.xori data.xori+data.xdom])
ylim([data.yori data.yori+data.ydom])

% plotting
myFolder = 'Figs';

width = data.sizecmf(1);   % cm
height = data.sizecmf(2);   % cm
set(gcf,'paperunits','centimeters')
set(gcf,'units','centimeters')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize',[width,height])
set(gcf,'paperposition',[0,0,width,height])
if ~exist(myFolder, 'dir')
    mkdir(myFolder)
end
% adding f to the name
myRealname = [data.myFigname,'s'];
myfile = fullfile(myFolder, myRealname);
%print(myfile,'-dpdf')
print(myfile,'-depsc')
%print(myfile,'-dpng')
hold off

end

% find new point in zero contour at a distance dl
function [xp,dphi,lambda]=findnewpoint(data,nmax,tolmin,xo,xp,lambda, dl, dvar)

% initialize
mytol = 1.;
it = 1;

while mytol>tolmin && it<nmax
    % B surface
    [phi, dphi, d2phi]=bspline2d(xp,data,dvar);
    
    % right hand side
    dxp = (xp - xo)';
    rhs = [- phi*dphi(1) + lambda*dxp(1);...
        - phi*dphi(2) + lambda*dxp(2);...
        - dxp(1)*dxp(1) - dxp(2)*dxp(2) + dl^2];
    %Amat = dphi'*dphi + phi * [d2phi(1) d2phi(2);d2phi(2) d2phi(3)] - lambda*eye(2)
    Kcol1 = [dphi(1)*dphi(1) + phi*d2phi(1) - lambda; ...
        dphi(1)*dphi(2) + phi*d2phi(2);...
        2*dxp(1)];
    Kcol2 = [ dphi(1)*dphi(2) + phi*d2phi(2); ...
        dphi(2)*dphi(2) + phi*d2phi(3) - lambda;...
        2*dxp(2)];
    Kcol3 = [ -dxp(1); ...
        -dxp(2);...
        0];
    % solution
    %deltax = Kmat\rhs;
    
    % Cramer's rule
    %detK = det([Kcol1 Kcol2 Kcol3])
    detK = Kcol1(2)*Kcol2(3)*Kcol3(1)...
        +Kcol1(3)*Kcol2(1)*Kcol3(2)...
        -Kcol1(3)*Kcol2(2)*Kcol3(1)...
        -Kcol1(1)*Kcol2(3)*Kcol3(2);
    
    %detx1 = det([Kcol Kcol2 Kcol3]);
    detx1 = rhs(2)*Kcol2(3)*Kcol3(1)...
        +rhs(3)*Kcol2(1)*Kcol3(2)...
        -rhs(3)*Kcol2(2)*Kcol3(1)...
        -rhs(1)*Kcol2(3)*Kcol3(2);
    %detx2 = det([Kcol1 Kcol Kcol3]);
    detx2 = Kcol1(2)*rhs(3)*Kcol3(1)...
        +Kcol1(3)*rhs(1)*Kcol3(2)...
        -Kcol1(3)*rhs(2)*Kcol3(1)...
        -Kcol1(1)*rhs(3)*Kcol3(2);
    %detx3 = det([Kcol1 Kcol2 Kcol]);
    detx3 = Kcol1(1)*Kcol2(2)*rhs(3)...
        +Kcol1(2)*Kcol2(3)*rhs(1)...
        +Kcol1(3)*Kcol2(1)*rhs(2)...
        -Kcol1(3)*Kcol2(2)*rhs(1)...
        -Kcol1(1)*Kcol2(3)*rhs(2)...
        -Kcol1(2)*Kcol2(1)*rhs(3);
    
    %     plot(xp(1),xp(2),'r*')
    %     hold on
    
    % update
    xp(1) = xp(1) + detx1/detK;
    xp(2) = xp(2) + detx2/detK;
    lambda = lambda + detx3/detK;
    mytol = abs(phi);
    it =it+1;
end
if mytol>tolmin
    'No convergence for new point'
    xp
end
end

function [xp,dphi] = startingp(xa,ya,xb,yb, data, tolmin, dvar)

% point a
a = 0.;
phia=bspline2d([xa,ya],data,dvar);

% point b
b = 1.;
phib=bspline2d([xb,yb],data,dvar);

x0 = xa;
y0 = ya;
dx = xb-xa;
dy = yb-ya;
tolmy = 1;
xp=[xa, ya];
dphi = [0,0];

nsame = 10;
da = 1/(nsame + 1);
nit= 1;
% bisect up to ntimes too if no sign change this is not
while (nit <= nsame && phia*phib>0) %limit iterations to prevent infinite loop and solution found
    % new first point
    a = b - nit*da;
    xa = x0 + a*dx;
    ya = y0 + a*dy;
    fprintf('No sign change: new point a %6.3f %6.3f \n',xa,ya)
    xp = [xa, ya];
    [phia,dphi,~] = bspline2d(xp,data,dvar);
    nit = nit + 1; %increment step counter
    
end


% check there is a sign change
if phia*phib <=0
    nit = 1;
    nmax = 100;
    while (nit <= nmax && tolmy > tolmin) %limit iterations to prevent infinite loop and solution found
        % new midpoint
        c = (a + b)/2;
        xc = x0 + c*dx;
        yc = y0 + c*dy;
        xp = [xc, yc];
        [phic,dphi,~] = bspline2d([xc,yc],data,dvar);
        tolmy = abs(phic);
        nit = nit + 1; %increment step counter
        % new interval
        if phia*phic > 0
            a = c;
            phia = phic;
        else
            b = c;
        end
    end
    if tolmy>tolmin
        'No convergence for starting point'
        %         thetao = thetao + pi/50;
        %         ntry = ntry + 1;
        %         %re initialize
        %         mytol = 1.;
        %         it = 1;
    end
else
    'No sign change for bisection alg to get starting point'
end
end

function [phi,dphi,d2phi,ivar,NB,dNBx,dNBy,dNBxx,dNBxy,dNByy]=bspline2d(xpts,data,rho)
xori=data.xori;
yori=data.yori;
xdom=data.xdom;
ydom=data.ydom;
nexv=data.nexv;
neyv=data.neyv;

% middle
xm = xori + .5*xdom;
ym = yori + .5*ydom;

xpt = xpts;

%size of elements
dxv=(xdom-xori)/nexv;
dyv=(ydom-yori)/neyv;

%segment b spline Determine the design variable number
ixv=max(min(ceil((xpt(1)-xori)/dxv),nexv),1);
iyv=max(min(ceil((xpt(2)-yori)/dyv),neyv),1);

%index of element node for the element (in first layer)

ivar=(neyv+3)*(ixv-1)+iyv-1+...
    [1, 2, 3, 4, ...
    neyv+4, neyv+5, neyv+6, neyv+7,...
    2*neyv+7, 2*neyv+8, 2*neyv+9,2*neyv+10,...
    3*neyv+10,3*neyv+11,3*neyv+12,3*neyv+13];

%center point of segment
xmid = 0.5d0*dxv+xori + (ixv-1)*dxv;
ymid = 0.5d0*dyv+yori + (iyv-1)*dyv;

%isoparametric mapping
r = 2.d0*(xpt(1)-xmid)/dxv;
s = 2.d0*(xpt(2)-ymid)/dyv;

%B-SPLINE shape functions
Neb=[r.^3 r.^2 r 1]*[-1 3 -3 1;3 -3 -3 3;-3 -15 15 3;1 23 23 1]/48;
Nnb=[s.^3 s.^2 s 1]*[-1 3 -3 1;3 -3 -3 3;-3 -15 15 3;1 23 23 1]/48;

%Shape functions in 2D
NB=[Nnb*Neb(1) Nnb*Neb(2) Nnb*Neb(3) Nnb*Neb(4)];

%scalar function
phi=NB*rho(ivar);

%derivatives
%B-SPLINE shape functions
dNeb=[3*r.^2 2*r 1 0]*[-1 3 -3 1;3 -3 -3 3;-3 -15 15 3;1 23 23 1]/48;
dNnb=[3*s.^2 2*s 1 0]*[-1 3 -3 1;3 -3 -3 3;-3 -15 15 3;1 23 23 1]/48;

%Shape functions in 2D
dNBx=[Nnb*dNeb(1) Nnb*dNeb(2) Nnb*dNeb(3) Nnb*dNeb(4)]*2/dxv;

%derivative respect to control points
dphidx=dNBx*rho(ivar);

%Shape functions in 2D
dNBy=[dNnb*Neb(1) dNnb*Neb(2) dNnb*Neb(3) dNnb*Neb(4)]*2/dyv;

%derivative respect to control points
dphidy=dNBy*rho(ivar);

dphi=[dphidx dphidy];

% second derivatives
%B-SPLINE shape functions
ddNeb=[6*r 2 0 0]*[-1 3 -3 1;3 -3 -3 3;-3 -15 15 3;1 23 23 1]/48;
ddNnb=[6*s 2 0 0]*[-1 3 -3 1;3 -3 -3 3;-3 -15 15 3;1 23 23 1]/48;


% %Second derivative Shape functions
% Ee=Se\(ddNB-Le*Be)
dNBxx = [Nnb*ddNeb(1) Nnb*ddNeb(2) Nnb*ddNeb(3) Nnb*ddNeb(4)]*4/(dxv^2);
dNBxy = [dNnb*dNeb(1) dNnb*dNeb(2) dNnb*dNeb(3) dNnb*dNeb(4)]*4/(dxv*dyv);
dNByy = [ddNnb*Neb(1) ddNnb*Neb(2) ddNnb*Neb(3) ddNnb*Neb(4)]*4/(dyv^2);

d2phidxx = dNBxx*rho(ivar);
d2phidxy = dNBxy*rho(ivar);
d2phidyy = dNByy*rho(ivar);
d2phi=[d2phidxx d2phidxy d2phidyy];

end

% plot Bsurface
function plotBsurf(data,dvarM)

% vector form of design variables
dvar = dvarM(:);

%unpack
xori=data.xori;
yori=data.yori;
xdom=data.xdom;
ydom=data.ydom;
Xcs =data.Xcs;
Ycs =data.Ycs;
nexv=data.nexv;
neyv=data.neyv;

%size of elements
dxv=(xdom-xori)/nexv;
dyv=(ydom-yori)/neyv;

%surface domain
nplot=100;
xplot=linspace(xori,xdom,nplot);
yplot=linspace(yori,ydom,nplot);

%mesh of domain
[X,Y]=meshgrid(xplot,yplot);
figure(1)
%empty Z
Z=zeros(size(X));

%bsurface
for i=1:nplot
    for j=1:nplot
        [phi]=bspline2d([X(i,j),Y(i,j)],data,dvar);
        
        %bsurface val
        Z(i,j)=phi;
    end
    
end
% levels
v = [0.,0.];
%surf of control points
[C,h]=contour(X,Y,Z,v,'LineWidth',1);
axis equal
hold on
xlabel('$x_1$')
ylabel('$x_2$')
text(.1,.9,'$\Omega$','FontSize',data.myFontSizeLabel)
text(.7,.45,'$\partial \Omega_0$','FontSize',data.myFontSizeLabel)

% store as a pdf
% plotting
sizecmf = [8 8];
myFolder = 'Figs';

width = sizecmf(1);   % cm
height =sizecmf(2);   % cm
set(gcf,'paperunits','centimeters')
set(gcf,'units','centimeters')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize',[width,height])
set(gcf,'paperposition',[0,0,width,height])
if ~exist(myFolder, 'dir')
    mkdir(myFolder)
end
% adding f to the name
myRealname = [data.myFigname,'f'];
myfile = fullfile(myFolder, myRealname);
%print(myfile,'-dpdf')
print(myfile,'-depsc')
%print(myfile,'-dpng')
hold off

nCpoi = C(2,1);
% center point and mean radius of hole
data.xc = sum(C(:,2:end)')/nCpoi;
data.rc = sum(sqrt(C(1,2:end).^2 + C(2,2:end).^2))/nCpoi;
figure()
h=surf(X, Y, Z);
h.EdgeColor = 'interp';
h.LineStyle = '-';
h.LineWidth = 0.1;
%make it just a mesh
shading interp
alpha 0.9
lightangle(-45,30)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.3;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
hold on
%surf of control points
h=surf(Xcs,Ycs,dvarM);


%code patch('Vertices',V,'Faces',F,'FaceColor',[1 0 0])
V = [xori-dxv yori-dyv 0;xori-dxv yori+ydom+dyv 0;xori+xdom+dxv yori+ydom+dyv 0;xori+xdom+dxv yori-dyv 0];
F = [1,2,3,4];
patch('Vertices',V,'Faces',F,'FaceColor',[1 0 0],'Facealpha',.5,'Edgealpha',0.5)
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\phi$')
set(h,'edgecolor','k','facecolor','none','linewidth',.1)
hold on
s=10; %100*ones(size(dvarM(:)));

%control points as circles
h=scatter3(Xcs(:),Ycs(:),dvarM(:),s,'MarkerEdgeColor',[0 .7 .7],...
    'MarkerFaceColor','k',...
    'LineWidth',1);
%axis equal
%,...
%   'MarkerEdgeColor','k',...
%  'MarkerFaceColor',[0 .75 .75],'MarkerSize',1)
width = sizecmf(1);   % cm
height =sizecmf(2);   % cm
set(gcf,'paperunits','centimeters')
set(gcf,'units','centimeters')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize',[width,height])
set(gcf,'paperposition',[0,0,width,height])
if ~exist(myFolder, 'dir')
    mkdir(myFolder)
end
% adding f to the name
myRealname = [data.myFigname,'c'];
myfile = fullfile(myFolder, myRealname);
%print(myfile,'-dpdf')
print(myfile,'-depsc')
%print(myfile,'-dpng')
end

