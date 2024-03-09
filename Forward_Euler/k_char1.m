function kchar1=k_char1(Ts)
%Limitador de temperatura
Ts = min(Ts, 610);
kchar1=2.66E10*exp(-12800/Ts); %Se calcula la constante de reaccion del modelo de formacion de
                                                            %carbonizado a partir de la biomasa presente en la cama de combustible [1/s].
                                                            %Esta constante de reaccion fue tomada de https://doi.org/10.1016/S0010-2180(99)00124-8