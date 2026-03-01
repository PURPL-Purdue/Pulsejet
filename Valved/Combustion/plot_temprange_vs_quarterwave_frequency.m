T_range = 500:10:2500;      %temp range from 500 to 2500K
gamma = 1.4;                %adiabatic index
R = 287;                    %specific gas constant air (J/kg*K)

c = sqrt(R*T_range*gamma);        %speed of sound at specific temperature (m/s)
L = 1.83;                         %6ft = 1.83m
f_range = c/(4*L);                %frequency calculation

plot(T_range, f_range)
xlabel("Temperature (K)")
ylabel("Frequency (Hz)")
title("Temp vs. Frequency for quarter-wave model")
grid()