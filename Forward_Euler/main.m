clear all; close all; clc;
%
L=4; 
H=0.2; 
b=2; 
vel_parrilla=0.0014;    %Velocidad de avance de la parrilla movil [m/s]. 

%Discretizacion temporal
t_f=2000;   %Tiempo final de la simualcion [s].
delta_t=10E-3; %Incremento de tiempo usado para el avance temporal [s].
deltat_save=10;  %Incremento de tiempo para guardado de datos [s].
n_deltat=round(t_f/delta_t);   %Numero de pasos de tiempo totales.
n_deltat_save=round(deltat_save/delta_t); %Numero de pasos de tiempo entre registros.
t=0;   %Valor del tiempo inicial [s]
tplot=0; %Vector de tiempos para ploteo

%Discretizacion espacial
n_nod=20; %Numero de divisiones en direccion z [-].
delta_z=H/(n_nod-1);  %Distancia entre nodos [m]. Alto del volumen de control.
xplot=linspace(0,H,n_nod);  %Vector de coordenadas nodales para ploteo
%delta_x=vel_parrilla*delta_t; %Ancho del volumen de control [m].
delta_x=0.5;
vol_cell=delta_x*delta_z*b;   %Volumen de cada celda [m3].

T_ing=293.15;              %Temperatura de ingreso del combustible al calcinador [K]. Tomado de http://dx.doi.org/10.1108/HFF-04-2018-0165
                                        %En http://dx.doi.org/10.1016/j.fuel.2017.06.037 tambien se reporta una temperatura ambiente del combustible de la cama a la entrada de la parrilla. 

T_furn=1200;               %Temperatura interna de los gases del horno [K]. Tomado de http://dx.doi.org/10.1016/j.fuel.2012.08.017
                                       %En http://dx.doi.org/10.1108/HFF-04-2018-0165 se plantea una temperatura de reactor de 1123 K (pag. 516).                                       
%Aire primario
p_atm= 74660.5;          %Presion atmosferica [Pa]
Tair_in=293.15;            %Temperatura del aire primario de alimentacion [K].
v_air=0.2;                    %Velocidad del aire primario de alimentacion [m/s].

%Biosolido
rho_lodo=1100;           %Densidad del lodo [kg/m3].
porc_moist=0.2;          %Porcentaje de humedad inicial de los lodos que ingresan al calcinador [-].
rho_s=densidad_solidos(rho_lodo,porc_moist); %Densidad del solido seco [kg/m3].
c_s=500;                      %Calor especifico de la fase solida [J/kg/K]. 
                                      %Valor de 2000 [J/kg/K] en http://dx.doi.org/10.1016/j.fuel.2017.06.037
k_s=0.2;                       %Coeficiente de conduccion de calor [W/m/K]. Tomado de http://dx.doi.org/10.1016/j.fuel.2017.06.037
dp=10E-3;                    %Diametro de la particula de solidos del lodo [m]. 
Nu=4.364;                    %Numero de Nusselt [-].

%Postproceso
z1=H;                            %Coordenada 1 para analisis
z2=0.5*H;                     %Coordenada 2 para analisis
z3=0.25*H;                   %Coordenada 3 para analisis

%Inicializacion de vectores
%Se incializan los vectores de almacenamiento
Ts=zeros(n_nod,1);
Tg=zeros(n_nod,1);
 M=zeros(n_nod,10,1);

%Se inicializa los vectores para la integracion temporal
Ts_temp_n=zeros(n_nod,1);
Tg_temp_n=zeros(n_nod,1);
Ts_temp_n1=zeros(n_nod,1);
Tg_temp_n1=zeros(n_nod,1);
dMdt_temp=zeros(n_nod,10);
M_temp_n1=zeros(n_nod,10);
dMdt=zeros(n_nod,10);
Qs=zeros(n_nod,1);
Qg=zeros(n_nod,1);

%Definicion de condiciones iniciales
Ts(:,1)=T_ing;
Tg(:,1)=T_ing;
Ts_temp_n(:,1)=T_ing;
Tg_temp_n(:,1)=T_ing;
pro=porc_moist*rho_lodo;
M(:,1,1)=pro;       %Concentracion inicial de humedad en la cama combustble [kg/m3].
M(:,2,1)=(1-porc_moist)*rho_lodo; %Concentracion inicial de biomasa en la cama de combustble [kg/m3].
M(:,3,1)=0;                                         %Concentracion inicial de volatiles en la cama de combustble.
M(:,4,1)=0;                                         %Concentracion inicial de carbonizado en la cama de combustble.
M(:,5,1)=0;                                         %Concentracion inicial de alquitranes en la cama de combustble.
M(:,6,1)=0;                                         %Concentracion inicial de CO en la cama de combustble.
M(:,7,1)=0;                                         %Concentracion inicial de CO2 en la cama de combustble.
M(:,8,1)=0;                                         %Concentracion inicial de CH4 en la cama de combustble.
M(:,9,1)=0;                                         %Concentracion inicial de O2 en la cama de combustble.
M(:,10,1)=0;                                       %Concentracion inicial de vapor de agus en la cama de combustble.
M_temp_n=M(:,:,1);

%Se definen las condiciones de borde para las especies en z=0
rho_aire=rho_g(Tair_in,p_atm);
m_O2=0.2*rho_aire*1000/32;  %Moles de O2 en el aire de entrada [mol/m3].
BC_g=[0 0 0 0 m_O2 0];

count=0;
entradas=1;
for i=1:n_deltat %Se itera en el tiempo
    %Se hace un barrido por cada nodo de la columna de combustible
    for j=1:n_nod
         [dMdt_temp,Qs(j),Qg(j)]=dM_dt(t,M_temp_n(j,:),Ts_temp_n(j),Tg_temp_n(j),pro,p_atm);
         dMdt(j,:)=dMdt_temp;
         M_temp_n1(j,1:4)=dMdt(j,1:4) * delta_t + M_temp_n(j,1:4);  %Se calcula el nuevo valor de masa de especies en la fase solida
    end
    Ts_temp_n1=compute_Ts(Ts_temp_n,Tg_temp_n,delta_z,delta_t,k_s, rho_s,c_s,Nu,dp,Qs,Tair_in,T_furn);
    Tg_temp_n1=compute_Tg(Tg_temp_n,Ts_temp_n,delta_z,delta_t,p_atm,Nu,dp,Qg,Tair_in,v_air);

    C_temp=compute_C(M_temp_n(:,5:10),delta_z,delta_t,BC_g,dp,v_air,dMdt(:,5:10));
    M_temp_n1(:,5:10)=C_temp;
    
    count=count+1;
     t=i*delta_t;
     if count==n_deltat_save
        tplot=[tplot t];
        msg=strcat('tiempo= ', num2str(tplot(end)),' [s]');
        disp(msg);
        %Se almacenan los vectores de temperatura 
        entradas=entradas+1;
        Ts(:,entradas)=Ts_temp_n1;
        Tg(:,entradas)=Tg_temp_n1;
        %Se almacena los datos de masa
        M(:,:,entradas)=M_temp_n1;
        %Se grafican las temperaturas de la fase solida y gaseosa
        subplot(2,2,1);
        plot(xplot,Ts(:,entradas),'-r','LineWidth',1.5);
        ylim([273 1300]); xlabel('z [m]'); ylabel('T [K]');
        hold on; grid on;
        plot(xplot,Tg(:,entradas),'-b','LineWidth',1.5); 
        legend('Fase solida','Fase gaseosa');
        drawnow;
        count=0;
        % pause
        hold off;
     end
    Ts_temp_n=Ts_temp_n1;
    Tg_temp_n=Tg_temp_n1;
    M_temp_n=M_temp_n1;
end
tplot=tplot/3600;  %Se convierte el vector de tiempos de [s] a [h]
%Se grafica la temperatura de la fase solida y gaseosa a diferentes alturas a lo largo del tiempo
nod1=round(z1/delta_z);
nod2=round(z2/delta_z);
nod3=round(z3/delta_z);
Ts_nod1=Ts(nod1,1:end);
Ts_nod2=Ts(nod2,1:end);
Tg_nod1=Tg(nod1,1:end);
Tg_nod2=Tg(nod2,1:end);
subplot(2,2,2);
plot(tplot,Ts_nod1,'-b','LineWidth',1.5); 
hold on; grid on;
plot(tplot,Ts_nod2,'--b','LineWidth',1.5);
plot(tplot,Tg_nod1,'-r','LineWidth',1.5);
plot(tplot,Tg_nod2,'--r','LineWidth',1.5);
xlabel('Tiempo [h]');
ylabel('Temperatura [K]');
msg1=strcat('T_s [K] - z=',num2str(delta_z*round(z1/delta_z,2)),' [m]');
msg2=strcat('T_s [K] - z=',num2str(delta_z*round(z2/delta_z,2)),' [m]');
msg3=strcat('T_g [K] - z=',num2str(delta_z*round(z1/delta_z,2)),' [m]');
msg4=strcat('T_g [K] - z=',num2str(delta_z*round(z2/delta_z,2)),' [m]');
leg=legend(msg1,msg2,msg3,msg4,'Location','southoutside');
leg.NumColumns = 2;

%Se grafica el compotamiento de las masas en el tiempo
M1_nod1=zeros(entradas,1);
M2_nod1=zeros(entradas,1);
M3_nod1=zeros(entradas,1);
M4_nod1=zeros(entradas,1);
M5_nod1=zeros(entradas,1);
M6_nod1=zeros(entradas,1);
M7_nod1=zeros(entradas,1);
M8_nod1=zeros(entradas,1);
M9_nod1=zeros(entradas,1);
M10_nod1=zeros(entradas,1);
M1_nod2=zeros(entradas,1);
M2_nod2=zeros(entradas,1);
M3_nod2=zeros(entradas,1);
M4_nod2=zeros(entradas,1);
M5_nod2=zeros(entradas,1);
M6_nod2=zeros(entradas,1);
M7_nod2=zeros(entradas,1);
M8_nod2=zeros(entradas,1);
M9_nod2=zeros(entradas,1);
M10_nod2=zeros(entradas,1);
for i=1:entradas
    M1_nod1(i)=M(nod1,1,i);
    M2_nod1(i)=M(nod1,2,i);
    M3_nod1(i)=M(nod1,3,i);
    M4_nod1(i)=M(nod1,4,i);
    M5_nod1(i)=M(nod1,5,i);
    M6_nod1(i)=M(nod1,6,i);
    M7_nod1(i)=M(nod1,7,i);
    M8_nod1(i)=M(nod1,8,i);
    M9_nod1(i)=M(nod1,9,i);
    M10_nod1(i)=M(nod1,10,i);
    M1_nod2(i)=M(nod2,1,i);
    M2_nod2(i)=M(nod2,2,i);
    M3_nod2(i)=M(nod2,3,i);
    M4_nod2(i)=M(nod2,4,i);
    M5_nod2(i)=M(nod2,5,i);
    M6_nod2(i)=M(nod2,6,i);
    M7_nod2(i)=M(nod2,7,i);
    M8_nod2(i)=M(nod2,8,i);
    M9_nod2(i)=M(nod2,9,i);
    M10_nod2(i)=M(nod2,10,i);
end
msg1=strcat(' z=',num2str(delta_z*round(z1/delta_z,2)),' [m]');
msg2=strcat(' z=',num2str(delta_z*round(z2/delta_z,2)),' [m]');
msg3=strcat(' z=',num2str(delta_z*round(z1/delta_z,2)),' [m]');
msg4=strcat(' z=',num2str(delta_z*round(z2/delta_z,2)),' [m]');

%Se grafican las masas M1 a M4 - Fase solida
subplot(2,2,3);
yyaxis left;
plot(tplot,M1_nod1,'-b','LineWidth',1.5); 
hold on; grid on; 
plot(tplot,M1_nod2,'--b','LineWidth',1.5); 
plot(tplot,M2_nod1,'-r','LineWidth',1.5); 
plot(tplot,M2_nod2,'--r','LineWidth',1.5);
plot(tplot,M4_nod1,'-k','LineWidth',1.5); 
plot(tplot,M4_nod2,'--k','LineWidth',1.5); 
yyaxis right;
plot(tplot,M3_nod1,'-g','LineWidth',1.5); 
plot(tplot,M3_nod2,'--g','LineWidth',1.5); 
xlabel('t [h]'); 
leg=legend(strcat('Humedad [kg/m^3]',' - ', msg1),strcat('Humedad [kg/m^3]',' - ', msg2),...
             strcat('Biosolido [kg/m^3]',' - ', msg1),strcat('Biosolido [kg/m^3]',' - ', msg2),...
             strcat('Carbonizado [kg/m^3]',' - ',msg1),strcat('Carbonizado [kg/m^3]',' - ',msg2),...
             strcat('Volatiles [kg/m^3]',' - ',msg1),strcat('Volatiles [kg/m^3]',' - ',msg2), ...
             'Location','southoutside');
leg.NumColumns = 2;

%Se grafican las masas M5 a M10 - Fase gaseosa
subplot(2,2,4);
plot(tplot,M5_nod1,'-b','LineWidth',1.5); 
hold on; grid on;
plot(tplot,M5_nod2,'--b','LineWidth',1.5); 
plot(tplot,M6_nod1,'-r','LineWidth',1.5); 
plot(tplot,M6_nod2,'--r','LineWidth',1.5);
plot(tplot,M7_nod1,'-k','LineWidth',1.5); 
plot(tplot,M7_nod2,'--k','LineWidth',1.5);
plot(tplot,M8_nod1,'-g','LineWidth',1.5); 
plot(tplot,M8_nod2,'--g','LineWidth',1.5);
xlabel('t [h]'); 
leg=legend(strcat('CH_x O_y [mol/m^3]',' - ', msg1),strcat('CH_x O_y [mol/m^3]',' - ', msg2),...
             strcat('CO [mol/m^3]',' - ', msg1),strcat('CO [mol/m^3]',' - ', msg2),...
             strcat('CO_2 [mol/m^3]',' - ',msg1),strcat('CO_2 [mol/m^3]',' - ',msg2),...
             strcat('CH_4 [mol/m^3]',' - ',msg1),strcat('CH_4 [mol/m^3]',' - ',msg2), ...
             'Location','southoutside');
leg.NumColumns = 2;
