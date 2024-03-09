%Esta funcion define las ecuaciones diferenciales empleadas para modelar 
%el proceso de secado de la biomasa y la tasa de transformacion de la biomasa 
%en volatiles y carbonizado durante el proceso de pirolisis. Adicionamente 
%considera los procesos de oxidacion de CO, CH4 y CHO (alquitranes). La base
%del modelo fue tomada del articulo: Numerical modeling of municipal waste bed incineration. 
%DOI: http://dx.doi.org/10.1108/HFF-04-2018-0165.
%Los parametros de entrada son:
% ->t: tiempo 
% ->Mn: Vector de concentracion de componentes y especies en el tiempo n: 
%         (1) humedad remanente en la cama de combustible [kg/m3].
%         (2) biomasa residual en la cama de combustible [kg/m3].
%         (3) volatiles formados a partir de la biomasa [kg/m3].
%         (4) carbonizado [kg/m3].
%         (5) concentracion de alquitranes por unidad de volumen (CHO) [mol/m3].
%         (6) concentracion de CO por unidad de volumen [mol/m3]. 
%         (7) concentracion de CO2 por unidad de volumen [mol/m3].
%         (8) concentracion de CH4 por unidad de volumen [mol/m3].
%         (9) concentracion de O2 por unidad de volumen [mol/m3].
%         (10) concentracion de H2O por unidad de volumen [mol/m3].
% ->Ts: Temperatura de la fase solida [K]
% ->Tg: Temperatura de la fase gaseosa [K]
% ->initial_moist: masa inicial de agua contenida en el lodo dentro del volumen de control [kg/m3].
% ->p_atm: presion atmosferica [Pa].
%Los parametros de salida son:
%<-dMdt: Vector con las tasas de formacion o consumo de las especies [mol/m3/s] o [kg/m3/s].
%<-Q_s: Calor total generado o consumido por las diferentes reacciones sobre la fase solida [W/m3].
%<-Q_g: Calor total generado o consumido por las diferentes reacciones sobre la fase gaseosa [W/m3].
%
%Autor: Carlos Galeano. Universidad Nacional de Colombia.
%
function [dMdt,Q_s,Q_g]=dM_dt(t,Mn,Ts,Tg,initial_moist,p_atm)
%% Parametros del proceso
alfa=0.2835;                %Fraccion masica de volatiles que corresponde a masa de gas. Tomado de https://doi.org/10.1016/j.fuel.2012.08.017.
beta=0.7165;               %Fraccion masica de volatiles que corresponde a masa de alquitran. Tomado de https://doi.org/10.1016/j.fuel.2012.08.017.
gamma_CO2=0.65;    %Fraccion masica de los gases que corresponde a la masa de CO2. Tomado de https://doi.org/10.1108/HFF-04-2018-0165.
gamma_CO=0.3;        %Fraccion masica de los gases que corresponde a la masa de CO. Tomado de https://doi.org/10.1108/HFF-04-2018-0165.
gamma_CH4=0.05;    %Fraccion masica de los gases que corresponde a la masa de CH4. Tomado de https://doi.org/10.1108/HFF-04-2018-0165.
m_C=12/1000;            %Masa molar del carbono [kg/mol].
m_O2=32/1000;         %Masa molar del oxigeno [kg/mol].
m_CO2=44/1000;      %Masa molar del dioxido de carbono [kg/mol].
m_CO=28/1000;         %Masa molar del monoxido de carbono [kg/mol].
m_CHO=95/1000;      %Masa molar del alquitran [kg/mol].
m_CH4=16.04/1000;  %Masa molar del metano [kg/mol]
m_H2O=18/1000;       %Masa molar del agua [kg/mol].
dMdt=zeros(1,10);
Q_s=0;
Q_g=0;
Mn2=Mn(2);

Mn4=Mn(4);
Mn5=Mn(5);
Mn6=Mn(6);

Mn8=Mn(8);
Mn9=Mn(9);
Mn10=Mn(10);
%% Reaccion de secado del lodo de desecho
kdry=k_dry(Ts); %Se calcula la constante de reaccion del modelo de secado.
secado=kdry*initial_moist*exp(-kdry*t);
dMdt(1)=-secado;  %Se define la ecuacion diferencial que modela el secado de la cama de combustible, [kg/m3/s].

Q_s = Q_s - abs(dMdt(1))*40.65E3/m_H2O; %El valor de 40.65 kJ/mol se tomo de: http://dx.doi.org/10.1108/HFF-04-2018-0165.                                                                           
%% Biomasa residual
% Consumo de biomasa
kvol=k_vol(Ts);
kchar1=k_char1(Ts);                        %Se calcula la constante de reaccion del modelo de formacion de carbonizado a partir de la biomasa presente en la cama de combustible [1/s].
dMdt(2)=-(kvol+kchar1)*Mn2;      %Variacion de biomasa en el cama de combustible. El primer termino 
                                                             %se relaciona con el consumo de biomasa para la formacion de volatiles. 
                                                             %El segundo termino se relaciona con el consumo de biomasa para la
                                                             %formacion de carbonizado, [kg/m3/s].
Q_s = Q_s + kchar1*Mn2*418E3;%El valor de 418 kJ/kg fue tomado de http://dx.doi.org/10.1108/HFF-04-2018-0165.
%Q_g = Q_g + kvol*Mn(2)*418E3;     %El valor de 418 kJ/kg fue tomado de http://dx.doi.org/10.1108/HFF-04-2018-0165.
%% Reacciones de pirolisis y combustion
% Formacion de volatiles
dMdt(3)=kvol*Mn2;                       %Tasa de formacion de volatiles a partir de biomasa [kg/m3/s]. 
%Una segunda alternativa para la reaccion de formacion de volatiles se plantea en 
%https://doi.org/10.1016/j.fuel.2004.09.020

%Formacion y consumo de carbonizado
%Se calculan los parametros de la reaccion de oxidacion del carbonizado
[coef_C,coef_O2,coef_CO,coef_CO2,molar_fraction_O2]=char_param_oxidat(Ts); 
 pO2=p_atm*molar_fraction_O2;
kchar2=k_char2(Ts,pO2); 
kporm4=kchar2*Mn4; %Se calcula la constante de reaccion del proceso de oxidacion del carbonizado [1/s]. 
dMdt(4)=kchar1*Mn2-kporm4;      %Formacion de carbonizado. El primer tiempo es formacion por 
                                                                         %pirolisis a partir de biomasa y el segundo consumo por combustion 
                                                                         %del mismo [kg/m3/s]. 
Q_s = Q_s + kporm4*111E3/m_C; %El valor de 111 kJ/mol fue tomado de http://dx.doi.org/10.1108/HFF-04-2018-0165 y 
                                                                        %coincide con el valor reportado en la Tabla 1 de https://doi.org/10.1021/ef950193p
%% Formacion y consumo de especies
%  Formacion y consumo de alquitranes (CH1.84O0.96)
k4=beta*kvol;                                                        %Constante de reaccion para la formacion de alquitranes a partir de la pirolisis de biomasa.
Te=0.5*(Ts+Tg);
k5=2.9E5*Te*exp(-9650/Te);                              %Constante de reaccion para el consumo de alquitranes debido a la combustion de alquitranes.    
combualqui=k5*Mn5^0.5*Mn9;                                                                                 %Este termino se tomo de https://doi.org/10.1021/ef950193p aunque tambien se emplea en 
                                                                                 %https://doi.org/10.1016/j.fuel.2004.09.020.
dMdt(5)=k4*Mn2/m_CHO...                              %Tasa de formacion de alquitranes. Es el 71.65% en masa de los volatiles formados.  
                 -combualqui;                          %El segundo termino corresponde al consumo del alquitran por la combustion del mismo.

%Q_g = Q_g - (k5*Mn(5)^0.5*Mn(9))*42E3;        %El valor de 42 kJ/mol fue tomado de http://dx.doi.org/10.1108/HFF-04-2018-0165.
% Formacion y consumo de CO
k6=kchar2*coef_CO/(m_C*coef_C);                %Coeficiente de reaccion basado en la fraccion masica CO/C que se presenta en la combustion del carbonizado.
                                                                                 %Por cada mol de C que se consume se produce 2*(1-1/theta) moles de CO.
k7=alfa*gamma_CO*kvol/m_CO;                      %Coeficiente de reaccion para la formacion de CO a partir de los gases de volatiles.
k8=3.25E7*exp(-15098/Tg);                                %Coeficiente de reaccion para el consumo de CO producto de la oxidacion de esta molecula.
oxicO=k8*Mn6*Mn9^0.5*Mn10^0.5; 
%Esta expresion fue tomada de https://doi.org/10.1016/j.fuel.2004.09.020.
k9=1.6E10*exp(-24157/Tg);                                %Coeficiente de reaccion para el consumo de CH4 producto de su combustion. Por cada mol de CH4 que se 
  combuch4=k9*Mn8^0.7*Mn9^0.8;                                                                                 %consume se produce una mol de CO.
                                                                                 %Esta expresion fue tomada de https://doi.org/10.1016/j.fuel.2004.09.020.
dMdt(6)=k6*Mn(4)...                                              %El primer termino es formacion de CO por combustion de carbonizado.
                +k7*Mn2...                                            %El segundo termino corresponde al aporte de CO proveniente de los gases de los volatiles.
                -oxicO...    %El tercer termino corresponde al consumo de CO en la reaccion de oxidacion de esta molecula.
                +combuch4...                 %El cuarto termino corresponde a la produccion de CO por la combustion de CH4.                                                     
                +combualqui;                           %El quinto termino corresponde a la produccion de CO por la combustion de alquitranes.                                                     
%Q_g = Q_g + (k8*Mn(6)*Mn(9)^0.5*Mn(10)^0.5)*284E3;   %El valor de 284 kJ/mol fue tomado de http://dx.doi.org/10.1108/HFF-04-2018-0165.
% Formacion de CO2
k10=kchar2*coef_CO2/(m_C*coef_C);               %Coeficiente de reaccion basado en la fraccion masica CO2/C que se presenta en la combustion del carbonizado.
k11=alfa*gamma_CO2*kvol/m_CO2;                  %Coeficiente de reaccion para la formacion de CO2 a partir de los gases de volatiles
dMdt(7)=k10*Mn(4)...                                             %Formacion de especie de CO2. El primer termino es formacion de CO2 por combustion de carbonizado.
                +k11*Mn2...                                           %El segundo termino corresponde al aporte de los gases de los volatiles. 
                +oxicO;       %El tercer termino corresponde a la formacion de CO2 en la reaccion de oxidacion de CO.
% Formacion y consumo de CH4
k12=alfa*gamma_CH4*kvol;                                  %Coeficiente de reaccion para la formacion de CH4 a partir de los gases de volatiles.
dMdt(8)= k12*Mn2/m_CH4...                               %El primer termino corresponde al aporte de CH4 proveniente de los gases de los volatiles.
                   -combuch4;                    %El segundo termino corresponde al consumo del CH4 por combustion del mismo.
%Q_g = Q_g + (k9*Mn(8)^0.7*Mn(9)^0.8)*520E3; %El termino 520 kJ/mol se tomo de http://dx.doi.org/10.1108/HFF-04-2018-0165.
% Consumo de O2
dMdt(9)= -1.5*(combuch4)...                  %El primer termino corresponde al consumo de O2 por la combustion de CH4. Por cada mol de CH4 que se oxida, se consume 1.5 moles de O2.
                   -0.5*(oxicO)...  %El segundo termino corresponde al consumo de O2 por la oxidacion de CO. Por cada mol de CO que se oxida, se consume 0.5 moles de O2.
                   -0.48*(combualqui);                       %El tercer termino corresponde al consumo de O2 por la oxidacion de alquitranes (CHO). Por cada mol de CHO que se oxida, 
                                                                                            %se consume 0.48 moles de O2.
% Consumo y formacion de vapor de agua
dMdt(10)=(secado)/m_H2O...     %El primer termino tiene que ver con el agua evaporada de la cama de combustible solido [mol/m3/s].                                                                                 
                   -0.5*oxicO...        %El segundo termino corresponde al consumo de agua en la oxidacion del CO
                    +2*(combuch4)...                     %El tercer termino se refiere a la formacion de agua en la combustion de CH4. Se forman dos moles de agua por cada mol de metano.
                    +0.92*(combualqui);                         %El cuarto termino se refiere a la formacion de agua en la combustion de alquitranes (CHO). Se forman 0.92 moles
                                                                                                %de agua por cada mol de alquitran.
