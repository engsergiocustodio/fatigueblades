%--------------------------------------------------------------------
%%%% "fatigue_turbine_m2.m" FUNCTION
%%%% IMPLEMENTED BY SERGIO CUSTODIO - engsergiocustodio@gmail.com
%--------------------------------------------------------------------

%MODELO 2 - Com turbulência
function [tensao,vida_ciclos,vida_anos] = fatigue_turbine_m2(tf,turbulence,vm,V_BETFEM,S_BETFEM,R,Su,C_wholer,m_wholer,vel_dinamic,omega_dinamic)
tensao=interp1(V_BETFEM,S_BETFEM,vm,'spline');

omega_vm=interp1(vel_dinamic,omega_dinamic,vm,'spline');


tempo=0:60/(2*omega_vm):tf;
sigma=ones(1,size(tempo,2))*tensao;
sigma=sigma.*cos(2*pi*omega_vm*tempo/60)*((1-R)/2)+((1+R)/2)*tensao;

sigma=sigma.*normrnd(1,turbulence,[1,size(tempo,2)]);
%%
rf2 = rainflow_turbine(sigma');
ciclos_rfturbine=sum(rf2(3,:));
sigma_f2=rf2(1,:)./(1-rf2(2,:)./Su);
Nwholer=((1/C_wholer)*sigma_f2).^(m_wholer);
Dano=rf2(3,:)./Nwholer;
Dano_total2=sum(Dano);
Dano_ano2=Dano_total2*365.25*24*60*60/tempo(end);

vida_ciclos=(omega_vm*tempo(end)/60)/Dano_total2;
vida_anos=1/Dano_ano2;


