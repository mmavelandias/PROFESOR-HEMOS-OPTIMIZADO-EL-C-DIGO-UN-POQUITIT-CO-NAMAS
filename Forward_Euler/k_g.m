%Esta funcion calcula la conductividad termica de la fase gaseosa de la cama de combustible
%en funcion de la temperatura. El parametro de entrada es:
%-> Tg: Temperatura del gas [K].
%El parametro de salida es:
%<- kg: Conductividad termica del gas [W/m/K]
%Formula tomada de la Tabla 4 de https://doi.org/10.1016/j.fuel.2004.09.020
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function kg=k_g(Tg)
kg=4.8E-4*Tg.^0.717; 