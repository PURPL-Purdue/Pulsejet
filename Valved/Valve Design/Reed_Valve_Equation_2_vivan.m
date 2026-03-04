clear; clc;
% --- Constants for Air ---
gamma = 1.4;
R = 287;                % Gas constant [J/kg*K]
p0 = 101325;            % Atmospheric pressure [Pa]
T0 = 293.15;            % Atmospheric temperature [K] (approx 20°C)
rho0 = p0 / (R * T0);   % Atmospheric density [kg/m^3]

% --- Parameters (Modify these based on your specific setup) ---
V = 0.05;               % Tank volume [m^3]
A = 0.0001;             % Orifice area [m^2]
tspan = [0 5];         % Time interval [s]
p_initial = 0.9 * p0;   % Initial tank pressure [Pa] (e.g., 10% of atm)

% --- Constants from the Paper ---
E = 200e9; % Modulus for steel (Pa)
w = 0.010; % Width of petal (10mm -> 0.01m)
l = 0.020; % Length of petal (20mm -> 0.02m)
I = (w * 0.0005^3) / 12; % Moment of inertia for 0.5mm thick petal


% Execute ode45
[t, p_out] = ode45(@(t, p) pdot(t, p), tspan, p_initial);


figure;
plot(t, p_out, '-o')
title('ode45() for pressure')
grid on;

function dpdt = pdot(t, p)
    % --- Constants for Air ---
    gamma = 1.4;
    R = 287;                
    p0 = 101325;            % Upstream (Atmospheric)
    T0 = 293.15;            
    rho0 = p0 / (R * T0);   
    
    % --- Fixed Valve/Tank Parameters (Define these inside or pass them in) ---
    V = 0.05;               
    A = 0.01;             
    l = 0.010;              % Valve length (m)
    w = 0.03;              % Valve width (m)
    E = 200e9;              % Modulus (Pa)
    I = 5e-12;              % Moment of inertia (m^4)
    b = 0.001;              % Overlap (m)
    d = 0.003;              % Hole diameter (m)

    % 1. Calculate the pressure ratio
    r = p / p0;

    % 2. Safety Check & Valve Deflection
    % Valve opens when Atmospheric (p0) > Tank (p)
    delta_p = p0 - p; 

    if delta_p > 0 && r < 1.0
        % Calculate Deflection
        q = delta_p * w; 
        y_bar = (q * l^4) / (8 * E * I);
        
        % Calculate Flow Coefficient (Loss)
        % Using max to prevent division by zero
        y_safe = max(y_bar, 1e-9);
        xi_sigma = 0.55 + 4*(b/d - 0.1) + 0.155 / ((y_safe/d)^2);
        
        % term0 represents the discharge coefficient (Cd) logic
        term0 = 1 / sqrt(1 + xi_sigma);
        assignin("base", "term0", term0)
    else
        % If p >= p0, flow stops and valve is closed
        dpdt = 0;
        return; 
    end

    % 3. Differential Equation Terms
    term1 = (gamma * p0) / (rho0^gamma);
    
    % Density factor inside the tank
    density_factor = ( (p / (T0 * R)) * (r^((1 - gamma) / gamma)) )^(gamma - 1);
    term2 = density_factor / V;

    % 4. Ideal Mass Flux (Theoretical)
    sqrt_inner = 2 * p0 * rho0 * (gamma / (gamma - 1)) * (r^(2/gamma)) * (1 - r^((gamma-1)/gamma));
    mdot_ideal = A * sqrt(max(0, sqrt_inner));

    % 5. Final dp/dt
    % Apply the flow coefficient (term0) to the mass flow
    dpdt = term1 * term2 * (term0 * mdot_ideal);
end

%, gamma, R, p1, T0, rho1, V, A