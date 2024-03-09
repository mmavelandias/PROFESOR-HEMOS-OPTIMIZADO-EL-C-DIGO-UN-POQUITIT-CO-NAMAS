%Esta funcion calcula el vector de temperatura de la fase sodila a lo largo de la columna de combustible 
%en el tiempo n+1.
%Los parametros de entrada son:
%->C_n: Matriz de concentracion de los gases en el tiempo n [mol/m3]. Las especies consideradas se arreglan por 
%             columnas y las filas de la matriz corresponden a cada nodo de la discretizacion. Las columnas son:
% (1) Alquitran (gas) CH1.84O0.96.
% (2) CO
% (3) CO2
% (4) CH4
% (5) O2
% (6) H2O
%->delta_z: Distancia entre nodos [m].
%->delta_t: Tamaño del intervalo de tiempo para la integracion temporal [s].
%->C_BC: Vector fila con los valores de las especies en z=0 (condicion de borde) [mol/m3].
%->dp: Diametro de particula del biosolido [m]. 
%->v_air: Velocidad del aire primario que entra por debajo de la parrilla movil [m/s]. 
%->r: Matriz de fuentes de masa asociadas a las transformaciones fisico-quimicas del combustible. 
%       Cada fila corresponde a un volumen de control y cada columna a una de las especies modeladas [mol/m3/s]. 
%El parametro de salida es 
%<-C: Matriz de concentraciones nodales en el tiempo n+1. Cada fila corresponde a un nodo y cada 
%columna a una de las especies modeladas [mol/m3]. 
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function C=compute_C(C_n,delta_z,delta_t,C_BC,dp,v_air,r)
n_nod=size(C_n,1);  %Numero de nodos de la discretizacion
%Se definen los valores del coeficiente de dispersion para cada gas 
%(difusion molecular) [m2/s]. Los valores son tomados de: 
%http://dx.doi.org/10.1016/j.fuel.2012.08.017.
Di=[18.2E-6, 18.07E-6, 13.81E-6, 19.52E-6, 18.2E-6, 21.78E-6];
Da=Di+0.5*v_air*dp; %Se calcula el coeficiente de dispersion para cada gas. La formula se 
                                     %tomo de: https://doi.org/10.1016/j.fuel.2004.09.020
C=zeros(n_nod,6); %Se incializa el vector de salida
beta=Da*delta_t/(delta_z^2);
gamma=v_air*delta_t/delta_z;
%Se hace un recorrido por cada nodo de la discretizacion (excepto el primero y el ultimo)
for i=2:n_nod-1
  for j=1:6
      C(i,j)=beta(j)*C_n(i+1,j)+(1-2*beta(j)-gamma)*C_n(i,j)+(beta(j)+gamma)*C_n(i-1,j)+r(i,j)*delta_t;
  end
end
%Se aplica la condicion de borde en el primer nodo (z=0).
C(1,:)=C_BC; 
%Se aplica la condicion de borde dC/dz=0 en el ultimo nodo (z=L).
C(n_nod,:)=C(n_nod-1,:);

