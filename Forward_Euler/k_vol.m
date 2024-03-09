function kvol=k_vol(Ts)
%Limitador de temperatura
Ts=min(Ts,610);
kvol=5.16E6*exp(-10700/Ts);       %Se calcula la constante de reaccion del modelo de formacion de
                                                           %volatiles a partir de la biomasa presente en la cama de combustible [1/s]. 
                                                           %Esta constante de reaccion fue tomada de https://doi.org/10.1016/S0010-2180(99)00124-8