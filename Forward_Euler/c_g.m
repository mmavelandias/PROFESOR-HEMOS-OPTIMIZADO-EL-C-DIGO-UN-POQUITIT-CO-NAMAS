%Esta funcion calcula el calor especifico de la fase gaseosa de la cama de combustible 
%en funcion de la temperatura. El parametro de entrada es:
%-> Tg: Temperatura del gas [K].
%El parametro de salida es:
%<- cg: Calor especifico del gas [J/kg/K]
%Formula tomada de la Tabla 4 de https://doi.org/10.1016/j.fuel.2004.09.020
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function cg=c_g(Tg)
cg=(0.99+1.22E-4*Tg-5.68E3*Tg.^-2)*1E3;   %[J/kg/K]
%cg=1.9327E-10*Tg.^4-7.9999E-7*Tg.^3+1.1407E-3*Tg.^2-4.489E-1*Tg+1.0575E3;  %[J/kg/K]
