clear all; close all; clc;
 
%INPUTS:

%velocity = 2.35;                         %fluid velocity
%omega=22.614*2*pi/60;                       %rad/s
%omega=25*2*pi/60;                       %rad/s
vel = [1.25:0.25:2.5];%3];
omega = [8.0 11.25 14.66 18.15 21.84 25.7];%1 29.73 33.92];
%omega = [8.0 11.25 13 16 18 20.5];%1 29.73 33.92];
omega = omega*2*pi/60;

vproj=[1.0 1.5 2.0 2.33 2.5];
omegaproj=interp1(vel,omega,vproj,'spline');
omegaproj(1)=0.001;
ITRMAX = 10000;
TOL    = 1e-4;
B=3;                              %Number of Blades
coord='NACA653618.txt';           %AIRFOIL COORDINATES [X,Y]
rho = 999;  %water 20ºC          %fluid density (Kg/m^3)
visc = 9.7937E-7; %water 20ºC    %fluid kinematic viscosity (m2/s)
reynolds = 3000000;
AOAs =20.0;   %degree              %Angle of Atack STALL
cd = 0.01;                    %inicial value
cl = 1.0;                    %inicial value

bladeshape = 'shape235.txt';
blade_shape = load(bladeshape);
polar = 'polar_3E6.txt';

% CALCULATE
r=blade_shape(:,1);
c=blade_shape(:,2);
b=blade_shape(:,3);


figure(1)
subplot(2,1,1)
h=plot(r,c,r,b); grid on;
set(h,'LineWidth',1.0);
%set(gca,'FontName','times','FontSize',16)
xlabel('Posição Radial - (m)')
%ylabel('Corda (m) e Ângulo de torção (rad)')
legend('Corda (m)','Ângulo de Torção (rad)','Location','Best');
%xlabel('Radial position - (m)')
%ylabel('Chrod and Twist')
%legend('Chord','Twist','Location','Best');
grid off
%figure(2)
subplot(2,1,2)
nprofiles=30;
npoints=15;
%plot_shape3Da(coord,npoints,nprofiles,r',c',b','elipses.txt');
plot_shape3D(coord,npoints,r',c',b');
xlabel('x - (m)')
ylabel('y - (m)')
zlabel('z - (m)')

%%
%TSR=[4:0.5:10];
%TSR=omega*r(end)./vel;
TSR=omegaproj*r(end)./vproj;
%vel=zeros(1,size(TSR,2));
CP=zeros(1,size(TSR,2));
POWER=zeros(1,size(TSR,2));
CE=zeros(1,size(TSR,2));
FN=zeros(size(r,1),size(TSR,2));
FT=zeros(size(r,1),size(TSR,2));
Re=zeros(size(r,1),size(TSR,2));
W=zeros(size(r,1),size(TSR,2));
%%
for i=1:size(TSR,2)
    i
%vel(i)=omega*r(end)/TSR(i);
%[POWER(i),CE(i),CP(i),W(:,i),FN(:,i),FT(:,i),MT(:,i)]=BET(B,bladeshape,polar,omega(i),vel(i),rho,reynolds,ITRMAX,TOL,AOAs);
%[FN(:,i),FT(:,i),MT(:,i)]=BET(B,bladeshape,polar,omega(i),vel(i),rho,reynolds,ITRMAX,TOL,AOAs);
[FN(:,i),FT(:,i),MT(:,i),CP(i),POWER(i)]=BET(B,bladeshape,polar,omegaproj(i),vproj(i),rho,reynolds,ITRMAX,TOL,AOAs);
%[FN(:,i),FT(:,i),MT(:,i),CP(i)]=BET(B,bladeshape,polar,omega(i),vel(i),rho,reynolds,ITRMAX,TOL,AOAs);
end

FR=((FN.^2+FT.^2).^(1/2));

  %% PLOT
   figure(3)
mesquita2014 = load('mesquita_2014.txt');   
%   %h = plot(wrev(vel),wrev(POWER)); grid on;
h = plot(vproj,POWER/1000,'-k',mesquita2014(:,1),mesquita2014(:,2),'xk'); grid on;

set(h,'LineWidth',1.0);
%set(gca,'FontName','times','FontSize',16)
legend('Este Trabalho','Mesquita et al. (2014) e Holanda et al. (2017)', 'LOCATION', 'BEST')
xlabel('Velocidade do Rio - (m/s)')
ylabel('Potência - (kW)')
xlim([1.0 2.5])
grid off
%    xlabel('TSR')
%    ylabel('CP')
%   
% %  for i=1:size(W,2)  
% %  Re(:,i) = W(:,i).*c/visc;
% %  end

 %%
 figure(4)
 for i=1:size(vproj,2)
 v=zeros(size(r,1),1);
 h=plot(r,FN(:,i)./1000);
 set(h,'LineWidth',1.0);
 hold on
 end
 legend('1,0 m/s','1,5 m/s','2,0 m/s','2,33 m/s','2,5 m/s','LOCATION','BEST')
xlabel('r - (m)');
ylabel('Força Normal - (kN/m)');

%%
 figure(5)
 for i=1:size(vproj,2)
 v=zeros(size(r,1),1);
 h=plot(r,FT(:,i)./1000);
 set(h,'LineWidth',1.0);
 hold on
 end
 legend('1,0 m/s','1,5 m/s','2,0 m/s','2,33 m/s','2,5 m/s','LOCATION','BEST')
xlabel('r - (m)');
ylabel('Força Tangencial - (kN/m)');

%%
 figure(6)
 for i=1:size(vproj,2)
 v=zeros(size(r,1),1);
 h=plot(r,MT(:,i)./1000);
 set(h,'LineWidth',1.0);
 hold on
 end
 legend('1,0 m/s','1,5 m/s','2,0 m/s','2,33 m/s','2,5 m/s','LOCATION','BEST')
xlabel('r - (m)');
ylabel('Momento Torçor - (kNm/m)');

%%
%fid = fopen('1_cargas0_5292.txt', 'w+');
%fid = fopen('1_5_cargas1_1781.txt', 'w+');
%fid = fopen('2_0_cargas1_9007.txt', 'w+');
%fid = fopen('2_33_cargas2_4150.txt', 'w+');
fid = fopen('2_5_cargas2_6923.txt', 'w+');

%for i=1:1%size(vel,2)
i=5;
for j=3:size(FR,1)
fprintf(fid,'%5.8f %5.8f %5.8f\n',FN(j,i)*(r(end)-r(end-1)),FT(j,i)*(r(end)-r(end-1)),MT(j,i)*(r(end)-r(end-1)));
end
%end 
fclose(fid);
 
