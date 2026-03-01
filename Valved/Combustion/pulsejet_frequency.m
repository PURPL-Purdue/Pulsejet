%frequency calculation 1/4 wave frequency model

T = input("Temperature (K): ");  %1000 K temperature assumed for now
gamma = 1.4;                %adiabatic index
R = 287;                    %specific gas constant air (J/kg*K)

c = sqrt(R*T*gamma);        %speed of sound at specific temperature (m/s)
L = input("Input length of pulsejet (m): ");     %6ft = 1.83m
f1 = c/(4*L);               %frequency calculation

fprintf("Frequency using 1/4 model: %.2f Hz \n", f1)

%using tan model (estimates between Hemholtz & 1/4 wave models)

V = input("Enter the ratio between the volume of combustion chamber and the tailpipe volume: "); %w/ dimensions from structures, V_chamber/V_tailpape = 1.314
w = pulsejet_omega(V);
Lt = .6*L;                   %tailpipe length, assuming .6 of total length (estimate)
f2 = (c/(2*pi*Lt)) * w;
fprintf("Frequency estimate between Helmholtz and 1/4 model: %.2f Hz \n", f2)

%chat helped with this:
function w = pulsejet_omega(Vbar)
    % Solves w*tan(w) = 1/Vbar for 0 < w < pi/2

    f = @(w) w.*tan(w) - 1./Vbar;

    % Bracket the physical root
    w_lower = 1e-6;
    w_upper = pi/2 - 1e-6;

    w = fzero(f, [w_lower, w_upper]);
end
