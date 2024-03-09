function kchar2=k_char2(Ts,pO2)
%Limitador de temperatura
Ts = min(Ts, 610);
               %Se calcula la presion parcial del oxigeno en la superficie de la
exponente=exp(-15900/Ts);                                                               %particula de carbonizado [Pa].
kchar2=8620*exponente*pO2;           %Se calcula la constante de reaccion del proceso de oxidacion del carbonizado [1/s]. 
                                                                       %Tomado de https://doi.org/10.1016/j.fuel.2004.09.020. 