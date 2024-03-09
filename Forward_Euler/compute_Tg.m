%Esta funcion calcula el vector de temperatura de la fase gaseosa a lo largo de la columna de combustible 
%en el tiempo n+1.
%Los parametros de entrada son:
%->Tg_n: Vector de temperatura de la fase gaseosa en la columna de combustible en el tiempo n [K].
%->Ts_n: Vector de temperatura de la fase solida en la columna de combustible en el tiempo n [K].
%->delta_z: Distancia entre nodos [m].
%->delta_t: Tamaño del intervalo de tiempo para la integracion temporal [s].
%->p_atm: Presion atmosferica [Pa] .
%->Nu: Numero de Nusselt [-].
%->dp: Diametro de las particulas que conforman la fase solida [m].
%->S: Relacion volumen a area de las particulas de la fase solida [m2/m3].
%->Q: Valor de la fuente de calor asociada a las transformaciones fisico-quimicas del combustible [W/m3]. 
%->T_air: Es la temperatura del aire primario que entra por debajo de la parrilla movil a la cama de combustible [K]. 
%->v_air: Velocidad del aire primario que entra por debajo de la parrilla movil [m/s]. 
%El parametro de salida es 
%<-Tg: Vector de temperaturas de la fase gaseosa sobre la columna de combustible calculadas para el tiempo n+1 [K].
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function Tg=compute_Tg(Tg_n,Ts_n,delta_z,delta_t,p_atm,Nu,dp,Q,T_air,v_air)
n_nod=size(Ts_n,1);  %Numero de nodos de la discretizacion
Tg=zeros(n_nod,1); %Se incializa el vector de salida
Tg(1)=T_air;  %Se impone la condicion de borde 
S=6/dp;  %Relacion volumen a area de las particulas de la fase solida [m2/m3].
%Se hace un recorrido por cada nodo de la discretizacion (excepto el primero y el ultimo)
for i=2:n_nod-1
    kg=k_g(Tg_n(i));                    %Se calcula el coeficiente de conduccion de calor a la temperatura del nodo [W/m/K]
    cg=c_g(Tg_n(i));                    %Se calcula el calor especifico a la temperatura del nodo [J/kg/K].
    rhog=rho_g(Tg_n(i),p_atm); %Se calcula la densidad del gas a la temperatura del nodo [kg/m3].
    gamma=delta_t/(rhog*cg);
    alfa=kg*gamma/(delta_z^2);
    delta=v_air*delta_t/delta_z;
    h=k_air(Tg_n(i))*Nu/dp;%Coeficiente convectivo entre la fase solida y la fase gaseosa en cada volumen de control [W/m2/K].
    beta=h*S*gamma;
   Tg(i)=alfa*Tg_n(i+1)+(1-delta-2*alfa-beta)*Tg_n(i)+(alfa+delta)*Tg_n(i-1)+beta*Ts_n(i)+Q(i)*gamma;
end
%Se aplica la condicion de borde para el ultimo nodo (z=L)
dT_dz=0;
Tg(end)=dT_dz*delta_z+Tg(end-1);
