%Esta funcion calcula la conductividad termica del aire en funcion de la
%temperatura. El parametro de entrada es:
%-> Tg: Temperatura del gas [K].
%El parametro de salida es:
%<- kg: Conductividad termica del aire [W/m/K]
%Formula tomada de la Tabla 4 de https://doi.org/10.1016/j.fuel.2004.09.020
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function kair=k_air(Tg)
kair=5.66E-5*Tg+1.1E-2; 