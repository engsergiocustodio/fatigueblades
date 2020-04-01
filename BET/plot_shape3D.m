%close all;  clear all;  format long
%----------------------------------
function [] = plot_shape3D(filename,npoints,radius,chord,twist)
XiYi = load(filename);
Ri = radius';
Ci = chord';
Bi = twist';

%----------------------------------------------------------------
    nx = length(XiYi(:,1));         nc = length(Ci);
    ny = length(XiYi(:,2));         nr = length(Ri);
%----------------------------------------------------------------
    minr = min(Ri);
    maxr = max(Ri);
%...............................................................
    [xm,ym] = baric_perfil(filename);

%...............................................................
%...............................................................
% INTERPOLA AS COORDENADAS DOS PONTOS DO PERFIL EM CADA ESTAÇAO
%---------------------------------------------------------------
   xmin = min(XiYi(:,1)-xm);       ymin = min(XiYi(:,2)-ym);
   xmax = max(XiYi(:,1)-xm);       ymax = max(XiYi(:,2)-ym);

   dx0 = (xmax-xmin)/(nx-1);       x0 = xmin:dx0:xmax;
   dy0 = (ymax-ymin)/(ny-1);       y0 = ymin:dy0:ymax;

   dx = (xmax-xmin)/(npoints-1);        x = xmin:dx:xmax;
   dy = (ymax-ymin)/(npoints-1);        y = ymin:dy:ymax;

   xx = spline(x0,XiYi(:,1)-xm,x);%%%%correction
   yy = spline(y0,XiYi(:,2)-ym,y);%%%%correction

   
%...............................................................
% ESTABELECE AS COORDENADAS DOS PONTOS DO PERFIL EM CADA ESTAÇAO
%---------------------------------------------------------------
             FX = Ci*xx;
			 FY = Ci*yy;

%.............................................................
% PROVOCA A ROTAÇAO DOS VETORES QUE FORMAM AS ESTAÇOES DA PA
%-------------------------------------------------------------
       
for i = 1:npoints
	        FXr(:,i) = cos(Bi).*FX(:,i) - sin(Bi).*FY(:,i);
            FYr(:,i) = sin(Bi).*FX(:,i) + cos(Bi).*FY(:,i);
end



z = linspace(minr,maxr,nr);

r = sqrt(FXr.*FXr + FYr.*FYr);

theta = atan2(FYr,FXr);

    Ax = r.*cos(theta);
    Ay = r.*sin(theta);
    Az = Ri;


mesh(Ax,Az,Ay); axis equal;
%set(gca,'FontName','times','FontSize',14)
xlabel('x-(m)');
ylabel('y-(m)');
zlabel('z-(m)');
axis on; grid off;  box on;
colormap([0 0 0]);

rotate3d;


function [xm,ym] = baric_perfil(filename)
%-------------------------------
XiYi = load(filename);
%-------------------------------
    nx = length(XiYi(:,1));
    ny = length(XiYi(:,2));
%-------------------------------

x = XiYi(:,1);
y = XiYi(:,2);

for i = 1:nx-1
    Ax(i)  = (x(i+1)-x(i))*y(i);
    Ay(i)  = (y(i+1)-y(i))*x(i);
    xA(i) = (x(i+1)+x(i))/2*Ax(i);
    yA(i) = (y(i+1)+y(i))/2*Ay(i);
end

xm = sum(xA)/sum(Ax);
ym = sum(yA)/sum(Ay);

