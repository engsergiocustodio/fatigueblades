%--------------------------------------------------------------------
%%%% "fadiga_maincode"
%%%% IMPLEMENTED BY SERGIO CUSTODIO - engsergiocustodio@gmail.com
%--------------------------------------------------------------------

clear all; close all; clc;
%% Tempo e Velocidade
tf=10000;  %tempo de analise
turbulence=0.00;%turbulencia em relacao a velocidade media
seriehistorica=load('seriehistorica.txt');
%seriehistorica=[1.30 1.0;2.33 0.00; 2.6 0.00];
%% Dinâmica

dinamica=load('dinamica.txt');
vel_dinamic = dinamica(:,1)';
omega_dinamic = dinamica(:,2)';

%% BET-FEM
BET_FEM=load('BET_FEM.txt');
Velocidade=BET_FEM(:,1)'; %m/s
Scasca=BET_FEM(:,2)';

%% S-N

Su=310; %MPa %tensao de ruptura %MPa %Aluminium AA6061-T6 MPa Zakaria 2013
% Curva de Wholer sigma_a_equi=C_wholer*(N_wholer1^(1/m_wholer))
C_wholer=650.8; %MPa %Aluminium AA6061-T6 MPa Zakaria 2013
m_wholer=-1/0.12;  %Aluminium AA6061-T6 Zakaria 2013

%% Resultados
R_lista = [0.5 0.33 0.1 -1.0];

%V_lista = [1.0 1.5 2.0 2.33 2.5]; 

Resultado = zeros(size(R_lista,2),8);

contador=1;
for k=1:size(seriehistorica,1)
for i=1:size(R_lista,2)
vm=seriehistorica(k,1);
R=R_lista(1,i);        %razao de tensoes
Resultado(contador,1)=vm;
Resultado(contador,4)=R;
[Resultado(contador,2),Resultado(contador,5),Resultado(contador,6)]=fatigue_turbine_m2(tf,turbulence,vm,Velocidade,Scasca,R,Su,C_wholer,m_wholer,vel_dinamic,omega_dinamic);
Resultado(contador,7)=seriehistorica(k,2);
Resultado(contador,3)=Resultado(contador,7).*0.01./Resultado(contador,6);
contador=contador+1;
end
end

%% Vida por Velocidade
figure(1)

for j=1:size(R_lista,2)
for i=1:size(seriehistorica,1)
cont=size(R_lista,2)*(i-1)+j;
vidaciclos(j,i)=Resultado(cont,5);
end
h=plot(seriehistorica(:,1),vidaciclos(j,:));
set(h,'LineWidth',1.0);
hold on
end
legend('R=0.5','R=0.33','R=0.1','R=-1.0','LOCATION','BEST')
set(h,'LineWidth',1.0);
set(gca, 'YScale', 'log')
set(gca,'ylim',[1e3 1e21])
set(h,'LineWidth',1.5);
xlabel('River Velocity (m/s)')
ylabel('Life (cicles)')
xlim([min(seriehistorica(:,1)) max(seriehistorica(:,1))])

figure(2)

for j=1:size(R_lista,2)
for i=1:size(seriehistorica,1)
cont=size(R_lista,2)*(i-1)+j;
vidaanos(j,i)=Resultado(cont,6);
end
h=plot(seriehistorica(:,1),vidaanos(j,:));
set(h,'LineWidth',1.0);
hold on
end
vidaanos(end+1,:)=30;
h=plot(seriehistorica(:,1),vidaanos(end,:),'--');
hold on
vidaanos(end+1,:)=10;
h=plot(seriehistorica(:,1),vidaanos(end,:),'--');
hold on
vidaanos(end+1,:)=1;
h=plot(seriehistorica(:,1),vidaanos(end,:),'--');
hold on

legend('R=0.5','R=0.33','R=0.1','R=-1.0','30 years','10 years','1 year','LOCATION','BEST')
set(h,'LineWidth',1.0);
set(gca, 'YScale', 'log')
set(gca,'ylim',[1e-5 1e15])
set(h,'LineWidth',1.5);
xlabel('River Velocity (m/s)')
ylabel('Life (years)')
xlim([min(seriehistorica(:,1)) max(seriehistorica(:,1))])
%%
figure(3)

for j=1:size(R_lista,2)
for i=1:size(seriehistorica,1)
cont=size(R_lista,2)*(i-1)+j;
danos(j,i)=Resultado(cont,3);
end
vidatotal(j,1)=1/sum(danos(j,:),2);
h=stem(seriehistorica(:,1),danos(j,:),'filled');
hold on
end

legend('R=0.5','R=0.33','R=0.1','R=-1.0','LOCATION','BEST')
set(gca, 'YScale', 'log')
xlabel('River Velocity (m/s)')
ylabel('Damage per year')
xlim([min(seriehistorica(:,1)) max(seriehistorica(:,1))])

vidatotal

