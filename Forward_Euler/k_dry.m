%Esta funcion permite calcular la constante de reaccion del modelo 
%de secado del biosolido.
%Expresion tomada del articulo: Numerical modeling of municipal waste 
%bed incineration. DOI: http://dx.doi.org/10.1108/HFF-04-2018-0165.
%El parametro de entrada es:
%->Ts: Temperatura de la fase solida de la cama de combustible [K].
%El parametro de salida es:
%<-kdry: Coeficiente de reaccion [1/s].
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function kdry=k_dry(Ts)
%Limitador de temperatura
Ts=min(Ts,371.15);

Ts=((Ts-273.15)*(9/5))+32;  %Se pasa la temperatura de K a F.
beta=1.14;
kdry=5.16*10^(-100*beta)*Ts^48;
%kdry=0.0158;