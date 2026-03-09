% This code calculates the adiabatic flame temperature (Max theoretical temp) inside the
% engine when using propane and air at different ratios. The ratio is represented with an
% equivalence ratio (actual fuel to air ratio / ideal fuel to air ratio), and the AFT is
% calculated by finding the moles of each product and reactant, and multiplying it by the
% specific enthalpy for each molecule, and setting enthalpy of reactants = enthalpy of products
clear;
clc;

% Constants
hf_C3H8 = -103847;
hf_O2   =       0;
hf_N2   =       0;
hf_CO2  = -393456;
hf_H2O  = -241845;
hf_CO   = -110541;

% NASA 7 high-temp coefficients (1000-5000K)
NASA.CO2 = [3.857460290E+00,  4.414370260E-03, -2.214814040E-06,  5.234901880E-10, -4.720841640E-14];
NASA.H2O = [2.677703880E+00,  2.973185670E-03, -7.737688890E-07,  9.443351400E-11, -4.268999750E-15];
NASA.N2  = [2.952576370E+00,  1.396900400E-03, -4.926316030E-07,  7.860101950E-11, -4.607552040E-15];
NASA.O2  = [3.612213960E+00,  7.485816030E-04, -1.982130790E-07,  3.374615740E-11, -2.390537170E-15];
NASA.CO  = [2.715185610E+00,  2.062527430E-03, -9.988257710E-07,  2.300530080E-10, -2.036477160E-14];

% Reference temperature
Tref = 298;

phi_values = 0.1:0.1:3;
T_ad = zeros(size(phi_values));

for i = 1:length(phi_values)
    phi = phi_values(i);
    x = 5/phi;

    [M_CO2, M_CO, M_O2] = getMoles(phi, x);

    H_react = (1 * (hf_C3H8 + 0)) ...
            + (x * (0 + 0)) ...
            + ((3.76 * x) * (0 + 0));

    H_prod = @(T) ...
        M_CO2  * (hf_CO2 + intCp(T, NASA.CO2) - intCp(Tref, NASA.CO2)) + ...
        3.76*x * (hf_N2  + intCp(T, NASA.N2)  - intCp(Tref, NASA.N2))  + ...
        4      * (hf_H2O + intCp(T, NASA.H2O) - intCp(Tref, NASA.H2O)) + ...
        M_CO   * (hf_CO  + intCp(T, NASA.CO)  - intCp(Tref, NASA.CO))  + ...
        M_O2   * (hf_O2  + intCp(T, NASA.O2)  - intCp(Tref, NASA.O2));

    T_ad(i) = fzero(@(T) H_prod(T) - H_react, 2100);
end

figure;
plot(phi_values, T_ad, 'LineWidth', 2);
grid on;
xlabel('\phi (Equivalence Ratio)');
ylabel('T_{ad} (K)');
title('Adiabatic Flame Temperature vs Equivalence Ratio (NASA coefficients)');

function result = intCp(T, a)
    R = 8.314;
    result = R * (a(1)*T + (a(2)/2)*T^2 + (a(3)/3)*T^3 + (a(4)/4)*T^4 + (a(5)/5)*T^5);
end

function [M_CO2, M_CO, M_O2] = getMoles(phi, x)
    if phi >= 1
        N = [1 1; 2 1];
        b = [3; 2*x - 4];
        sol = N\b;
        M_CO2 = sol(1);
        M_CO  = sol(2);
        M_O2  = 0;
    else
        M_CO2 = 3;
        M_CO  = 0;
        M_O2  = x - 5;
    end
end
