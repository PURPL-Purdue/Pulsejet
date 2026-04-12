function results = reed_valve_area()
% """
% Simulates the dynamic opening of a reed valve under a time-varying
% pressure difference using a lumped model based on the reference paper.
% Computes the reed tip opening x(t) and converts it into effective orifice
% area A_0(t), while enforcing one-sided opening.
%
% Author: Joaquin Alarcon
% """

clear; clc; close all;

% Parameters from the model structure
params.reedFaceArea = 338.4e-6; % effective pressure-loaded area of the reed [m^2]
params.dampingCoeff = 8.0; % damping term that resists reed motion [N*s/m] 
params.stiffness = 600.0; % effective reed stiffness against opening [N/m]  
params.effectiveMass = 0.518e-3; % lumped mass used in the reed dynamic model [kg]
params.deflectionFactor = 1.549; % correction factor relating pressure load to tip motion
params.initialOpening = 0.0; % Initial reed tip opening at t = 0 [m]
params.maxOpening = 1.9e-3; % Maximum allowed opening before mechanical stop [m] 
params.areaGain = 0.019; % converts tip opening x(t) into effective area [m]

% Time settings
timeSpan = [0, 0.03]; % simulation time interval [s]
initialState = [params.initialOpening; 0.0]; % [x; xdot] initial tip opening and tip velocity 

% Prescribed pressure difference input
params.deltaP = @(t) pressureDifferenceInput(t); % deltaP(t) used to drive reed motion

% Solve ODE
options = odeset('RelTol',1e-8,'AbsTol',1e-10); % solver tolerances for stable time integration
[t, state] = ode45(@(t, y) reedValveODE(t, y, params), timeSpan, initialState, options);

tipOpening = max(state(:,1), 0); % enforce non negative opening. the reed cannot penetrate through the valve seat
tipVelocity = state(:,2); % tip velocity history from the ODE solution

% Postprocess
n = numel(t); % number of recorded time points
deltaPHistory = zeros(n,1); % stores DeltaP(t) at each time point
loadMultiplierHistory = zeros(n,1); % stores pressures load multiplier X(x)
orificeAreaHistory = zeros(n,1); % stores effective orifice area A_o(t)

for i = 1:n
    deltaPHistory(i) = params.deltaP(t(i)); % pressure difference at current time
    loadMultiplierHistory(i) = pressureLoadMultiplier(max(tipOpening(i), 0)); % load correction factor from the paper
    orificeAreaHistory(i) = params.areaGain * max(tipOpening(i), 0); % A_o(t) = areaGain * x(t)
end

% Package results
results.time = t; % time vector returned by ODE45 [s]
results.tipOpening = tipOpening; % Non-negative reed tip opening history [m]
results.tipVelocity = tipVelocity; % Reed tip velocity history [m/s]
results.deltaP = deltaPHistory; % applied pressure-difference history [Pa]
results.loadMultiplier = loadMultiplierHistory; % Pressure load multiplier X(x)
results.orificeArea = orificeAreaHistory; % effective orifice area history [m^2]
results.params = params; % model parameters used in this run

assignin('base', 'reedAreaResults', results); % export the structure to the MATLAB base workspace for easy inspection.

% Plots
figure(1);
tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

subplot(3,1,1);
plot(t, deltaPHistory/1000, 'LineWidth', 2.2, 'Color', [0.00, 0.20, 0.45]);
grid on;
box on;
xlabel("Time (s)");
ylabel('\DeltaP (kPa)');
title('Example Pressure Difference');

subplot(3,1,2);
plot(t, tipOpening*1000, 'LineWidth', 2.2, 'Color', [0.55, 0.10, 0.10]);
grid on;
box on;
xlabel("Time (s)");
ylabel('x(t) (mm)');
title('Reed Tip Opening')

subplot(3,1,3);
plot(t, tipVelocity*1000, 'LineWidth', 2.2, 'Color', [0.45, 0.20, 0.55]);
grid on;
xlabel("Time (s)");
ylabel("v(t) (mm/s)")
title("Reed Tip Velocity")

figure(2);
plot(t, orificeAreaHistory*1e6, 'LineWidth', 2.2, 'Color', [0.10, 0.45, 0.25])
grid on
xlabel('Time (s)', 'FontSize', 13)
ylabel('A_o(t) (mm^2)', 'FontSize', 13)
title('Effective Orifice Area')

figure(3);
plot(deltaPHistory/1000, orificeAreaHistory*1e6, 'LineWidth', 2.2, 'color', [0.85, 0.50, 0.10]);
grid on;
ylabel("A_o(t) (mm^2)");
xlabel("\DeltaP (kPa)");
title("Dynamic Area-Pressure Response");
figure(4);
plot(tipVelocity*1000, tipOpening*1000,  'LineWidth', 2.2, 'Color', [0.15, 0.55, 0.60]);
grid on;
ylabel("Reed Tip Opening [mm]");
xlabel("Reed Tip Velocity [mm/s]");
title("Phase Portrait of Reed Tip Motion")

end

%% Function reedValveODE
function dState = reedValveODE(t, state, params)
% REEDVALVEODE
% Computes the time derivative of the reed-valve state for the ODE solver.
%
% Inputs:
%  t: Current time [s]
%  state: Current state vector:
%         state(1) = tip opening x [m]
%         state(2) = tip velocity xdot [m/s]
%  params: structure containing the model parameters and the prescribed
%          pressure-difference input DeltaP(t)
%
% Outputs:
%  dState: Time derivative of the state vector:
%          dState(1) = xdot [m/s]
%          dState(2) = xddot [m/s^2]
%
% Model:
%  Uses the lumped reed dynamics from the reference paper:
%  m * xddot + x * xdot + c1 * k * (x - x0) = c1 * DeltaP * A * X

x = state(1); % tip opening x [m]
xDot = state(2); % tip velocity xdot [m/s]

deltaP = params.deltaP(t); % evaluate the prescribed pressure difference at the current time
X = pressureLoadMultiplier(max(x, 0)); % compute the pressure load multiplied X(x) from the paper correlation

% Paper-based lumped reed equation. rearranged to solve for xddot.
xDDot = ( ...
    params.deflectionFactor * deltaP * params.reedFaceArea * X ...
    - params.dampingCoeff * xDot ...
    - params.deflectionFactor * params.stiffness * (x - params.initialOpening) ...
    ) / params.effectiveMass;

% Closed-seat limit: do not allow penetration below zero opening
if x <= 0 && xDDot <= 0
    dState = [0; 0];
    return
end
% Stopper limit: do not allow motion beyond max opening
if x >= params.maxOpening && xDDot >= 0
    dState = [0; 0];
    return
end
dState = [xDot; xDDot]; % standard state-space output for the ODE solver
end

%% function: pressureLoadMultiplier
function X = pressureLoadMultiplier(x)
% PRESSURELOADMULTIPLIER
% Computes the pressure load multiplier X(x) using the empirical
% correlation given in Eq.(14) of the reference paper.
%
% Input:
%  x: Reed tip opening [m]
%
% Output:
%  X: Dimensionless pressure load multiplier    

numerator = -55.64*x + 0.07711; % Numerator of the empirical correlation
denominator = x^3 + 0.7044*x^2 - 49.5*x + 0.08965; % Denominator of the empirical correlation

% Enforce a non-negative multiplier.
X = max(0, numerator / denominator);
end

%% function pressureDifferenceInput -- PLACEHOLDER CHANGE AS NEEDED
function deltaP = pressureDifferenceInput(t)
% Demo pressure signal with both positive and negative pressure difference.
% Positive DeltaP opens the reed; negative DeltaP helps close it.

meanPressure = 10e3; % [Pa]
amp1 = 45e3; % [Pa]
amp2 = 18e3; % [Pa]
freq1 = 80; % [Hz]
freq2 = 160; % [Hz]

% 
deltaP = meanPressure + amp1 * sin(2*pi*freq1*t) + amp2 * sin(2*pi*freq2*t + 0.7);
end