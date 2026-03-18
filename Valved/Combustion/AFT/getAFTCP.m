function T = getAFTCP(phi)
% getAFTCP  Adiabatic flame temperature for propane/air at a given equivalence ratio.
%   T = getAFTCP(phi)  returns the AFT in Kelvin for the given phi.
%   Dissociation is accounted for via equilibrium constants from NASA polynomials.

% --- Formation enthalpies (J/mol) ---
hf.C3H8 = -103847;
hf.CO2  = -393456;
hf.H2O  = -241845;
hf.CO   = -110541;
hf.O2   =       0;
hf.N2   =       0;
hf.H2   =       0;
hf.OH   =   38987;
hf.O    =  249170;
hf.H    =  218000;

% --- NASA 7 high-temp coefficients (1000-5000K) [a1-a7] ---
NASA.CO2 = [3.857460290E+00,  4.414370260E-03, -2.214814040E-06,  5.234901880E-10, -4.720841640E-14, -4.875916600E+04,  2.271638060E+00];
NASA.H2O = [2.677703880E+00,  2.973185670E-03, -7.737688890E-07,  9.443351400E-11, -4.268999750E-15, -2.988589920E+04,  6.882555710E+00];
NASA.N2  = [2.952576370E+00,  1.396900400E-03, -4.926316030E-07,  7.860101950E-11, -4.607552040E-15, -9.239486880E+02,  5.871887620E+00];
NASA.O2  = [3.612213960E+00,  7.485816030E-04, -1.982130790E-07,  3.374615740E-11, -2.390537170E-15, -1.197815100E+03,  3.670331500E+00];
NASA.CO  = [2.715185610E+00,  2.062527430E-03, -9.988257710E-07,  2.300530080E-10, -2.036477160E-14, -1.415187500E+04,  7.818687720E+00];
NASA.H2  = [2.991423490E+00,  7.000644440E-04, -5.633828890E-08, -9.231578780E-12,  1.582752060E-15, -8.350340190E+02, -1.355110190E+00];
NASA.OH  = [2.867472880E+00,  1.056504480E-03, -2.590827580E-07,  3.052186740E-11, -1.331958760E-15,  3.718857740E+03,  5.701640730E+00];
NASA.O   = [2.569420780E+00, -8.597411370E-05,  4.194845890E-08, -1.001777990E-11,  1.228336910E-15,  2.921757910E+04,  4.784338640E+00];
NASA.H   = [2.500000010E+00, -2.308429730E-11,  1.615619480E-14, -4.735152350E-18,  4.981973570E-22,  2.547365990E+04, -4.466829140E-01];

Tref = 298;
P    = 1;

% Bounds: [n_CO2, n_CO, n_H2O, n_H2, n_O2, n_OH, n_O, T]
lb = [0,  0,  0,  0,  0,  0,  0,  1000];
ub = [10, 10, 10, 10, 10, 10, 10, 4000];

vars0   = [3, 0.01, 4, 0.01, 0.01, 0.01, 0.01, 2100];
options = optimoptions('lsqnonlin', 'Display', 'off');

x   = 5 / phi;
sol = lsqnonlin(@(v) equilibrium(v, x, NASA, Tref, hf, P), vars0, lb, ub, options);
T   = sol(8);

end

% =========================================================================

function F = equilibrium(vars, x, NASA, Tref, hf, P)
    n_CO2 = vars(1); n_CO  = vars(2); n_H2O = vars(3);
    n_H2  = vars(4); n_O2  = vars(5); n_OH  = vars(6);
    n_O   = vars(7); T     = vars(8);

    n_N2  = 3.76 * x;
    n_tot = n_CO2 + n_CO + n_H2O + n_H2 + n_O2 + n_OH + n_O + n_N2;

    F(1) = n_CO2 + n_CO - 3;
    F(2) = 2*n_H2O + 2*n_H2 + n_OH - 8;
    F(3) = 2*n_CO2 + n_CO + n_H2O + n_OH + 2*n_O2 + n_O - 2*x;

    Kp   = getKp(T, NASA);
    F(4) = (n_CO/n_tot)*P * ((n_O2/n_tot)*P)^0.5 - Kp(1) * (n_CO2/n_tot)*P;
    F(5) = (n_H2/n_tot)*P * ((n_O2/n_tot)*P)^0.5 - Kp(2) * (n_H2O/n_tot)*P;
    F(6) = (n_OH/n_tot)*P * ((n_H2/n_tot)*P)^0.5 - Kp(3) * (n_H2O/n_tot)*P;
    F(7) = ((n_O/n_tot)*P)^2                      - Kp(4) * (n_O2/n_tot)*P;

    F(8) = energyBalance(n_CO2, n_CO, n_H2O, n_H2, n_O2, n_OH, n_O, n_N2, T, NASA, Tref, hf) / 1e6;
end

function F = energyBalance(n_CO2, n_CO, n_H2O, n_H2, n_O2, n_OH, n_O, n_N2, T, NASA, Tref, hf)
    H_react = hf.C3H8;
    H_prod  = n_CO2*(hf.CO2 + intCp(T,NASA.CO2) - intCp(Tref,NASA.CO2)) + ...
              n_N2 *(hf.N2  + intCp(T,NASA.N2)  - intCp(Tref,NASA.N2))  + ...
              n_H2O*(hf.H2O + intCp(T,NASA.H2O) - intCp(Tref,NASA.H2O)) + ...
              n_CO *(hf.CO  + intCp(T,NASA.CO)  - intCp(Tref,NASA.CO))  + ...
              n_O2 *(hf.O2  + intCp(T,NASA.O2)  - intCp(Tref,NASA.O2))  + ...
              n_H2 *(hf.H2  + intCp(T,NASA.H2)  - intCp(Tref,NASA.H2))  + ...
              n_OH *(hf.OH  + intCp(T,NASA.OH)  - intCp(Tref,NASA.OH))  + ...
              n_O  *(hf.O   + intCp(T,NASA.O)   - intCp(Tref,NASA.O));
    F = H_prod - H_react;
end

function Kp = getKp(T, NASA)
    Kp(1) = exp(-(gibbsRT(T,NASA.CO)  + 0.5*gibbsRT(T,NASA.O2) - gibbsRT(T,NASA.CO2)));
    Kp(2) = exp(-(gibbsRT(T,NASA.H2)  + 0.5*gibbsRT(T,NASA.O2) - gibbsRT(T,NASA.H2O)));
    Kp(3) = exp(-(gibbsRT(T,NASA.OH)  + 0.5*gibbsRT(T,NASA.H2) - gibbsRT(T,NASA.H2O)));
    Kp(4) = exp(-(2*gibbsRT(T,NASA.O) -     gibbsRT(T,NASA.O2)));
end

function G = gibbsRT(T, a)
    G = a(1)*(1-log(T)) - (a(2)/2)*T - (a(3)/6)*T^2 - (a(4)/12)*T^3 - (a(5)/20)*T^4 + a(6)/T - a(7);
end

function result = intCp(T, a)
    R      = 8.314;
    result = R*(a(1)*T + (a(2)/2)*T^2 + (a(3)/3)*T^3 + (a(4)/4)*T^4 + (a(5)/5)*T^5);
end