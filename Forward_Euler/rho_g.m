%Esta funcion calcula la densidad del aire en funcion de la temperatura. 
% El parametro de entrada es:
%-> Tg: Temperatura del gas [K].
%-> Patm: Presion atmosferica [Pa].
%El parametro de salida es:
%<- rhog: Densidad del aire [kg/m3].
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function rhog=rho_g(Tg,p_atm)
R=287;   %[J/kg/K].  
rhog=p_atm./(R*Tg);  %[kg/m3]
