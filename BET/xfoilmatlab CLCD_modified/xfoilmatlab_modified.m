%%%Interacao XFOIL-MATLAB - SERGIO CUSTODIO 2018

clear all; close all; clc

coord='NACA653618.txt'; %%%%Arquivo com as coordenadas do perfil
iter= 10000;               %%%%Numero máximo de interacoes no xfoil
Re= [3e6];   %%%Reynolds input do usuario
Mach=0;      %%%Mach
alphamin=-5.0; %%%Menor angulo de ataque
alphamax=90;  %%%Maior angulo de ataque
dalpha=5;     %%%Intervalo

%% inputs
alpha=(alphamin:dalpha:alphamax);
N_Re = size(Re,2)
% Re= linspace(3e4,1e5,size(alpha,2));
cl=zeros(N_Re,size(alpha,2));
cd=zeros(N_Re,size(alpha,2));
cm=zeros(N_Re,size(alpha,2));

perfil=importdata(coord);
cp=zeros(size(perfil,1),size(alpha,2));
cs=zeros(size(perfil,1),size(alpha,2));

for p=1:N_Re
p
for i=1:size(alpha,2)
   % Write xfoil command file
i   
fid = fopen('xfoil.inp','w+');
if (fid<=0),
  error([mfilename ':io'],'Unable to create xfoil.inp file');
end;

fprintf(fid,'%s\n','plop');
fprintf(fid,'%s\n','g');
fprintf(fid,'\n');
fprintf(fid,'%s %s\n','load ',coord);
fprintf(fid,'\n');
fprintf(fid,'%s\n','oper');
fprintf(fid,'%s %f\n','iter ',iter);
fprintf(fid,'%s %f\n','re ',Re(p));
fprintf(fid,'%s %f\n','mach ',Mach);
fprintf(fid,'%s\n','visc');
fprintf(fid,'%s\n','pacc');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%s %f\n','alfa ',alpha(i));
fprintf(fid,'%s\n','cpwr cpwr.txt');
fprintf(fid,'%s\n','dump dump.txt');
fprintf(fid,'%s\n','pwrt');
fprintf(fid,'%s\n','pwrt.txt');
fprintf(fid,'%s\n','y');
fprintf(fid,'\n');
fprintf(fid,'%s\n','quit');
fclose(fid);

 % execute xfoil
   [status,result] = system('xfoil.exe < xfoil.inp > xfoil.out');
  if (status~=0),
    disp(result);
  end;
    
    cpwr=importdata('cpwr.txt');  %% é coeficiente de pressão (distribuição de forças ao redor do perfil)
    dump=importdata('dump.txt');  %% faz parte da distribuição de pressão. 
    
    fid=fopen('pwrt.txt','r');
    while true
        tline=fgetl(fid);
        if ~isempty(strfind(tline,'alpha'));break;end
    end
    
    tline=fgetl(fid);
    tline=fgetl(fid);
    
    if ~ischar(tline)
        status=-1;
        cli=0/0;
        cdi=0/0;
        cmi=0/0;
    
    else
        status=0;
        coefficients= str2num(tline);
        cli= coefficients(2);
        cdi= coefficients(3);
        cmi= coefficients(5);
    end
    fclose(fid);
    
    cl(p,i)=cli;
    cd(p,i)=cdi;
    cm(p,i)=cmi;
    
end

end



for p=1:N_Re
for i=2:size(cl,2)-1
if isnan(cd(i))==1
    cd(p,i)=(cd(p,i+1)+cd(p,i-1))/2;
    cl(p,i)=(cl(p,i+1)+cl(p,i-1))/2;
    cm(p,i)=(cm(p,i+1)+cm(p,i-1))/2;
end
end

Reynolds = [alpha; cl; cd; cm]';
end
%%
interp_values = zeros(size(alpha,2),size(Reynolds,2));
interp_values(:,1) = alpha; 

for i=2:size(Reynolds,2)
   interp_values(:,i) = interp1(alpha,Reynolds(:,i),alpha)
end 
%%
fid = fopen('polar.txt', 'w+');
fprintf(fid, '%i %i\n',size(Re,2),size(alpha,2));
for i=1:size(Re,2)
fprintf(fid, '%f\n',Re(i));
for j=1:size(interp_values,1)
fprintf(fid,'%5.8f %5.8f %5.8f %5.8f\n',interp_values(j,1),interp_values(j,i+1),interp_values(j,size(Re,2)+i+1),interp_values(j,2*size(Re,2)+i+1));
end
end 
fclose(fid);

%%%%%%%
% %%
% figure(1)
% subplot(4,1,1)
% plot(alpha,cl,'k');title('Cl x alpha');grid on;
% subplot(4,1,2)
% plot(alpha,cd,'k');title('Cd x alpha');grid on;
% subplot(4,1,3)
% plot(alpha,cl./cd,'k');title('Cl/Cd x alpha');grid on;
% subplot(4,1,4)
% plot(alpha,cm,'k');title('Cm x alpha');grid on;





