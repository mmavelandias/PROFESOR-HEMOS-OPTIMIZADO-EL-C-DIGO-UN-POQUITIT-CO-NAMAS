%Esta funcion calcula la densidad de los solidos contenidos en el lodo que
%alimenta el calcinador. Los parametros de entrada son:
%->rho_lodo: Densidad del lodo [kg/m3].
%->porc_moist: Porcentaje de agua en el lodo [-].
%El parametro de salida es:
%<-rho_s: Densidad de los solidos en el lodo [kg/m3].
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function rho_s=densidad_solidos(rho_lodo,porc_moist)

rho_agua=1000;  %Densidad del agua [kg/m3].
rho_rel_lodo=rho_lodo/rho_agua; %Densidad relativa del lodo [-]
rho_rel_s=(1-porc_moist)/(1/rho_rel_lodo-porc_moist);  %Se calcula la densidad relativa del solido [-]
rho_s=rho_rel_s*rho_agua;  %Se calcula la densidad del solido [kg/m3]