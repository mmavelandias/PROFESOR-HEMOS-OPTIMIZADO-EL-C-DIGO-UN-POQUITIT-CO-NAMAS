%Esta funcion calcula los parametros de la reaccion de oxidacion del carbonizado
% producido por la pirolisis de la biomasa. Las expresiones fueron tomadas del articulo: 
%Numerical modeling of municipal waste bed incineration. 
%DOI: http://dx.doi.org/10.1108/HFF-04-2018-0165.
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function [coef_C,coef_O2,coef_CO,coef_CO2,molar_fraction_O2]=char_param_oxidat(Ts)
rc=12*exp(-3300/Ts); % Relacion de formacion CO/CO2 del mecanismo de combustion
inv_rc=1/rc;
theta=(1+inv_rc)/(0.5+inv_rc); %Calculo de la relacion estequiometrica
inv_theta=1/theta;
coef_C=1;
coef_O2=inv_theta;
coef_CO=2*(1-inv_theta);
coef_CO2=2*inv_theta-1;
moles_totales=coef_C+coef_O2+coef_CO+coef_CO2;
molar_fraction_O2=coef_O2/moles_totales;