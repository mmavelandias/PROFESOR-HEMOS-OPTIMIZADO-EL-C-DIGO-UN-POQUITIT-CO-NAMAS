%Esta funcion calcula el vector de temperatura de la fase sodila a lo largo de la columna de combustible 
%en el tiempo n+1.
%Los parametros de entrada son:
%->Ts_n: Vector de temperatura de la fase solida en la columna de combustible en el tiempo n [K].
%->Tg_n: Vector de temperatura de la fase gaseosa en la columna de combustible en el tiempo n [K].
%->delta_z: Distancia entre nodos [m].
%->delta_t: Tamaño del intervalo de tiempo para la integracion temporal [s].
%->k_s: Conductividad termica de la fase solida [W/m/K].
%->rho_s: Densidad de la fase solida (biomasa) de la cama de combustible [kg/m3].
%->c_s: Calor especifico de la fase solida (biomasa) de la cama de combustible [J/kg/K]
%->Nu: Numero de Nusselt [-].  
%->dp: Diametro de las particulas que conforman la fase solida [m].
%->Q: Vector de fuentes de calor asociadas a las transformaciones fisico-quimicas del combustible. 
%         Cada fila corresponde a un volumen de control [W/m3]. 
%->T_air: Es la temperatura del aire primario que entra por debajo de la parrilla movil a la cama de combustible [K]. 
%->T_furn:Es la temperatura al interior del calcinador [K].
%El parametro de salida es 
%<-Ts: Vector de temperaturas de la fsae solida sobre la columna de combustible calculadas para el tiempo n+1 [K].
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function Ts=compute_Ts(Ts_n,Tg_n,delta_z,delta_t,k_s, rho_s,c_s,Nu,dp,Q,T_air,T_furn)
sigma=5.67E-8; %Constante de Stephan-Boltzman [W/m2/K4].
e_s=0.85;   %Factor de emisividad de la superficie [-]. Tomado de http://dx.doi.org/10.1108/HFF-04-2018-0165
n_nod=size(Ts_n,1);  %Numero de nodos de la discretizacion
Ts=zeros(n_nod,1); %Se incializa el vector de salida
Ts(1)=T_air;  %Se impone la condicion de borde 
S=6/dp;  %Relacion volumen a area de las particulas de la fase solida [m2/m3].
%Se hace un recorrido por cada nodo de la discretizacion (excepto el primero y el ultimo)
for i=2:n_nod-1
    k=k_air(Tg_n(i));
    gamma=delta_t/(rho_s*c_s);
    h=k*Nu/dp;  %Coeficiente convectivo entre la fase solida y la fase gaseosa en cada volumen de control [W/m2/K].
    beta=h*S*gamma;
    Keff=k_s+4*sigma*Ts_n(i)^3;  %El valor de 4 se toma de http://dx.doi.org/10.1108/HFF-04-2018-0165
    alfa=Keff*gamma/(delta_z^2);
    Ts(i)=alfa*Ts_n(i+1)+(1-2*alfa-beta)*Ts_n(i)+alfa*Ts_n(i-1)+beta*Tg_n(i)+Q(i)*gamma;
end
%Se aplica la condicion de borde para el ultimo nodo (z=L)
Keff=k_s+4*sigma*Ts_n(end)^3;
k=k_air(Tg_n(end));
h=k*Nu/dp;  
dT_dz=(h*(T_furn-Ts_n(end))+sigma*e_s*(T_furn^4-Ts_n(end)^4))/Keff;
Ts(end)=dT_dz*delta_z+Ts(end-1);





