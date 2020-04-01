%--------------------------------------------------------------------
%%%% "fatigue_turbine_m1.m" FUNCTION
%%%% IMPLEMENTED BY SERGIO CUSTODIO - engsergiocustodio@gmail.com
%--------------------------------------------------------------------

%MODELO 1 - Sem turbulência

function [tensao,Nwholer,vida_anos] = fatigue_turbine_m1(vm,V_BETFEM,S_BETFEM,Fator_sigma_max,R,Su,C_wholer,m_wholer,vel_dinamic,omega_dinamic)

tensao=interp1(V_BETFEM,S_BETFEM,vm,'spline');

omega_vm=interp1(vel_dinamic,omega_dinamic,vm,'spline');

sigma_max=tensao*Fator_sigma_max;
sigma_min=R*sigma_max;
sigma_a=(sigma_max-sigma_min)/2;
sigma_m=(sigma_max+sigma_min)/2;
sigma_a_equi=sigma_a/(1-sigma_m/Su);
Nwholer=((1/C_wholer)*sigma_a_equi)^(m_wholer);
ciclosAno=omega_vm*60*24*365.25;
vida_anos=Nwholer/ciclosAno;
end

