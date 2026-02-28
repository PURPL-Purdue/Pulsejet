function simulate_pressure()
    % --- Constants for Air ---
    gamma = 1.4;
    R = 287;                % Gas constant [J/kg*K]
    p0 = 101325;            % Atmospheric pressure [Pa]
    T0 = 293.15;            % Atmospheric temperature [K] (approx 20°C)
    rho0 = p0 / (R * T0);   % Atmospheric density [kg/m^3]
    
    % --- Parameters (Modify these based on your specific setup) ---
    V = 0.05;               % Tank volume [m^3]
    A = 0.0001;             % Orifice area [m^2]
    tspan = [0 1.5];         % Time interval [s]
    p_initial = 0.9 * p0;   % Initial tank pressure [Pa] (e.g., 20% of atm)

    % --- ODE Solver ---
    % 'p' is the state variable representing tank pressure
    [t, p_out] = ode45(@(t, p) pressure_ode(t, p, gamma, R, p0, T0, rho0, V, A), tspan, p_initial);

    results = [t, p_out];

    % This 'pushes' the variable 'PressureData' to the Base Workspace
    assignin('base', 'PressureData', results)

    % --- Plotting Results ---
    figure;
    subplot(2,1,1)
    plot(t, p_out / 1000, 'LineWidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('Pressure (kPa)');
    title('Tank Pressure Charging over Time');


% --- PART 2: Calculate and Plot Mass Flow (m_dot) ---
    % Note: We use p_out directly here instead of fetching from Workspace
    ratio = p_out ./ p0;
    sqrt_term = 2 .* p0 .* rho0 .* (gamma ./ (gamma-1)) .* (ratio .^ (2 ./ gamma)) .* (1 - ratio .^ ((gamma - 1)./gamma));
    
    % Clean up any imaginary numbers from numerical noise near ratio = 1
    sqrt_term(sqrt_term < 0) = 0; 
    m_dot = A .* sqrt(sqrt_term);

    subplot(2,1,2)
    plot(t, m_dot, 'r-', 'LineWidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('Mass Flow (kg/s)');
    title('Mass Flow Rate through Orifice');


end

function dpdt = pressure_ode(~, p, gamma, R, p1, T1, rho1, V, A)
    % 1. Calculate the pressure ratio (downstream/upstream)
    r = p / p1;

    % 2. First term of the equation: (gamma * p1 / rho1^gamma)
    term1 = (gamma * p1) / (rho1^gamma);

    % 3. Middle term: ( (p / (T1*R)) * (p/p1)^((1-gamma)/gamma) )^(gamma-1) / V
    % Note: This represents the (density)^(gamma-1) factor inside the tank
    density_factor = ( (p / (T1 * R)) * (r^((1 - gamma) / gamma)) )^(gamma - 1);
    term2 = density_factor / V;

    % 4. Square root term: Mass flux through orifice (Subsonic flow)
    % A * sqrt(2 * p1 * rho1 * (gamma/(gamma-1)) * r^(2/gamma) * (1 - r^((gamma-1)/gamma)))
    sqrt_inner = 2 * p1 * rho1 * (gamma / (gamma - 1)) * (r^(2 / gamma)) * (1 - r^((gamma - 1) / gamma));
    
    % Numerical guard: if pressure ratio exceeds 1, flow stops
    if sqrt_inner <= 0
        mdot = 0;
    else
        mdot = A * sqrt(sqrt_inner);
    end

    % Resulting dp/dt
    dpdt = term1 * term2 * mdot;
end

