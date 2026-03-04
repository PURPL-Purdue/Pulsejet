function Reed_valve_validation()
    % Gas Properties
    gamma = 1.4; % Ratio of specific heats (Cp/Cv)
    R = 287.0; % Specific gas constant for air [J/(kg*K)]
    
    % Upstream (ambient) conditions
    P_up = 101325; % Upstream (ambient) pressure [Pa] 
    T_up = 293.15; % upstream (ambient) pressure [K]
    rho_up = P_up / (R*T_up); % upsteam density. Assuming it is an ideal gas [kg/m^3]
    
    % Tank and Orifice Geometry
    tankVolume = 0.05; % Volume of the tank. For a fixed inflow mdot, a smaller volume increase pressure faster [m^3]
    A_orifice  = 1e-4; % Effective open area of the orifice [m^2]
    
    % Isothermal tank temperature. 
    % (It is reasonable just for basic "sanity check" because it avoids modeling heat transfer)
    T_tank = T_up; % Tank temperature stays constant over time [K]
    
    % Simulation settings
    % The ODE solver will pick its own internal time steps adaptively
    tSpan = [0, 1.5]; % defining the time interval for the simulation [s]
    % We start the tank at 90% of the upsteam pressure. It will be
    % "slightly under-pressurized".
    P0_tank = 0.90 * P_up; % Initial pressure inside the tank [Pa]
    
    % Define the ODE function handle
    % The ODE function returns dP_tank/dt at each time step based on
    % current tank pressure, flow model and 
    % isothermal tank relation: dp/dt = (R*T_tank/V_tank) * mdot_in
    odeFun = @(t, P_tank) tankPressureODE_isothermal(t, P_tank, gamma, R, P_up, T_up, T_tank, tankVolume, A_orifice);
    
    % Solving the ODE using ode45
    % This will give use the time vector that was chosen by solver and the
    % pressure solved at eac time value t.
    [t, P_tank] = ode45(odeFun, tSpan, P0_tank);
    
    % Compute mass flow rate mdot at each recorded time point
    % This uses the SUBSONIC compressible orifice model (THERE IS NO
    % CHOCKING YET). 
    % Upstream is fixed
    mdot = arrayfun(@(P) orificeMassFlow_subsonic(P_up, T_up, P, A_orifice, gamma, R), P_tank); % Compute mass flow rate
    
    % Package results and exports to workspace
    results = [t, P_tank, mdot]; % Combining time, pressure, and mass flow rate into one matrix.
    assignin('base', 'TankChargingData', results);
    
     %% --- Plot results ---
        figure(1);
        subplot(2,1,1);
        plot(t, P_tank/1000, 'LineWidth', 2);
        grid on;
        xlabel('Time (s)');
        ylabel('Tank Pressure (kPa)');
        title('Tank Pressure (Isothermal)');
        subplot(2,1,2);
        plot(t, mdot, 'LineWidth', 2);
        grid on;
        xlabel('Time (s)');
        ylabel('Mass Flow Rate (kg/s)');
        title('Mass Flow Rate Through Orifice (Subsonic Model)');
end

function dPdt = tankPressureODE_isothermal(~, P_tank, gamma, R, P_up, T_up, T_tank, V_tank, A_orifice)
    % If tank pressure is at/above upstream, no inflow. (I'm not truly
    % sure whether we should include backwards flow).
    if P_tank >= P_up
        dPdt = 0;
        return;
    end
    % Calling the other function in order to get the mass flow rate.
    mdot_in = orificeMassFlow_subsonic(P_up, T_up, P_tank, A_orifice, gamma, R); 
    % Isothermal ideal gas in fixed volume: 
    % dp/dt = (R*T/V)*mdot
    % For an ideal gas in vessel of constant volume V and constant
    % Temperature T.
    dPdt = (R * T_tank / V_tank) * mdot_in;
end


function mdot = orificeMassFlow_subsonic(P_up, T_up, P_down, A, gamma, R)
    % Returns subsonic compressible mass flow rate through an orifice.
    % Uses ideal gas: rho = P/(R*T)
    if A <= 0 || P_down >= P_up
        mdot = 0;
        return;
    end
    rho_up = P_up / (R * T_up); 
    r = P_down / P_up;  % pressure ratio (down/up), 0<r<1 for inflow
    % Subsonic compressible orifice mass-flux term (written with rho_up)
    sqrtInner = 2 * P_up * rho_up * (gamma/(gamma - 1)) * (r^(2/gamma)) * (1 - r^((gamma - 1)/gamma));
    % Numerical guard (near r ~ 1 can cause tiny negative due to roundoff)
    if sqrtInner <= 0
        mdot = 0;
    else
        mdot = A * sqrt(sqrtInner);
    end
end

