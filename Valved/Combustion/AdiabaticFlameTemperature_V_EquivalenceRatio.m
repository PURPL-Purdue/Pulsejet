clear;
clc;

% constant directly from the HW
hf_C3H8 = -103847;
hf_O2   =       0;    Cp_O2  = 35.593;
hf_N2   =       0;    Cp_N2  = 33.707;
hf_CO2  = -393456;    Cp_CO2 = 56.205;
hf_H2O  = -241845;    Cp_H2O = 43.874;
hf_CO   = -110541;    Cp_CO  = 34.148;

% original temp
Tref = 298;
phi_values = 0.1 : 0.1 : 3;

% empty vector for final temp
T_ad = zeros(size(phi_values));

for i = 1:length(phi_values)

    phi = phi_values(i);
    x = 1/(phi * 0.042 * 4.76);
    
    if phi >= 1
        N = [1 1; 2 1];
        M = [3; 2*x - 4];
    
        sol = N\M;
        A = sol(1);
        B = sol(2);
        C = 0;
    end
    
    if phi < 1
        A = 3;
        B = 0;
        C = x - 2 - A;
    end
    
    H_react = (1 * (hf_C3H8 + 0)) ...
            + (5 * (0 + 0)) ...
            + (18.8 * (0 + 0));
    
    H_prod = @(T_ad) (A * (hf_CO2 + Cp_CO2 * (T_ad - Tref))) ...
                   + (3.76 * x * (0 + Cp_N2 * (T_ad - Tref))) ...
                   + (4 * (hf_H2O + Cp_H2O * (T_ad - Tref))) ...
                   + (B * (hf_CO + Cp_CO * (T_ad - Tref))) ...
                   + (C * (+ Cp_O2 * (T_ad - Tref)));
    % find when enthalpy for reaction and product is the same
    T_ad(i) = fzero(@(T) H_prod(T) - H_react, 2100);

end

figure;
plot(phi_values, T_ad, 'LineWidth', 2);
grid on;
xlabel('\phi (Equivalence Ratio)');
ylabel('T_{ad} (K)');
title('Adiabatic Flame Temperature vs Equivalence Ratio');