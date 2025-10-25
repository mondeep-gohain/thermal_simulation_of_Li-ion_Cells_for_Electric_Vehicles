clear all
close all
%charge-single cell
%% Thermal and Geometric Parameters
L = 0.1;         % Length of the cylindrical cell (m)
R = 0.01;        % Radius of the cylindrical cell (m)
k = 3;           % Thermal conductivity (W/mK)
rho = 2700;      % Density (kg/m^3)
c = 800;         % Specific heat (J/kg.K)
h = 10;          % Convective heat transfer coefficient (W/m^2K)
T_ambient = 25;  % Ambient temperature (°C)
C_rate = 5; 

%% Battery & CCCV Charging Parameters
% Battery model parameters (Shepherd-model–inspired)
E0 = 4.3;          % Constant voltage term (V)
K = 0.001;         % Polarization constant (V/Ah)
Q_max = 10;         % Battery capacity (Ah)
R_internal = 0.01; % Internal resistance (Ohm)

V_max = 4.2;       % Maximum charging voltage (V)

delta_S = -5;   %Entropy change (J/mol.K)
F = 96485;       %Faradays const.(C/mol)

% Charging current parameters
I_CC = Q_max * C_rate;          % Constant current during CC phase (A)
I_min = Q_max * (C_rate/10);       % Minimum current threshold during CV phase (A)

% Time constant for exponential decay in CV phase (in seconds)
tau_CV = 50;       % Adjust this value based on desired decay rate

% Define SOC threshold for switching to CV mode (realistic value: 80–90% SOC)
SOC_threshold = 0.85;   % Transition from CC to CV at 85% SoC

% Initial state-of-charge (SoC) (range 0 to 1)
SoC = 0.2;
Q_current = SoC * Q_max;  % Current capacity (Ah)

% Calculate initial battery voltage using the Shepherd-like model:
V_batt = E0 - K * ((Q_max - Q_current) / Q_current) * I_CC - R_internal * I_CC;
if V_batt > V_max
    V_batt = V_max;
end

%% Derived Thermal Simulation Parameters
Volume = pi * R^2 * L;   % Cell volume (m^3)

% Spatial grid parameters
Nr = 50;                % Number of grid points in r-direction
Nz = 100;               % Number of grid points in z-direction
Dr = R / (Nr - 1);      % Radial step size
Dz = L / (Nz - 1);      % Axial step size

% Time stepping parameters
Dt = 0.01;              % Time step (s)
total_time = 3000;       % Total simulation time (s)
Nt = total_time / Dt;   % Number of time steps

%% Initialize Temperature Field
T = T_ambient * ones(Nr, Nz); % Temperature field (°C)
T_new = T;

%% Preallocate Arrays for Tracking Variables
time_array = (0:Nt-1) * Dt;
current_array = zeros(1, Nt);
voltage_array = zeros(1, Nt);
Q_joule_array = zeros(1, Nt);
Q_entropy_array = zeros(1, Nt);
Q_gen_array = zeros(1, Nt);
T_center = zeros(1, Nt);
RHS_term_array = zeros(1, Nt);
SoC_array = zeros(1, Nt);

found_t1 = false; % Flag to capture the CC-to-CV transition time
t1 = NaN;
T_t1 = [];      % To store temperature field at t1
T_final = [];   % To store final temperature field
center_j = ceil(Nz/2);  % Middle index in z-direction

%% Time-Stepping Loop for CCCV Charging & Thermal Simulation
for n = 1:Nt
    current_time = n * Dt;  % Current simulation time (s)
    
    % Decide charging mode based on SoC
    if SoC < SOC_threshold
        % CC Phase: Constant Current charging
        charge_mode = 'CC';
        I = I_CC;
        % Update SoC (convert Dt from seconds to hours)
        SoC = SoC + (I * Dt / 3600) / Q_max;
        if SoC > 1, SoC = 1; end
        Q_current = SoC * Q_max;
        % Update battery voltage using the Shepherd-like model:
        V_batt = E0 - K * ((Q_max - Q_current) / Q_current) * I - R_internal * I;
        if V_batt > V_max
            V_batt = V_max;
        end
    else
        % CV Phase: Constant Voltage charging
        charge_mode = 'CV';
        V_batt = V_max;  % Hold voltage constant in CV phase
        % Exponential decay of current, not falling below I_min:
        I = max(I * exp(-Dt/tau_CV), I_min);
        % Update SoC in CV phase as well:
        SoC = SoC + (I * Dt / 3600) / Q_max;
        if SoC > 1, SoC = 1; end
    end
    
    % Record electrical parameters and SoC
    current_array(n) = I;
    voltage_array(n) = V_batt;
    SoC_array(n) = SoC;
    
    % Calculate heat generation (volumetric) due to Joule and entropy effects
    T_avg = mean(T(:));
    Q_joule_volumetric = (I^2 * R_internal) / Volume;
    % Using delta_S = -50 J/mol.K and Faraday's constant F = 96485 C/mol:
    Q_entropy_volumetric_avg = (I * T_avg * (delta_S) / F) / Volume;
    
    
    Q_joule_array(n) = Q_joule_volumetric * Volume;
    Q_entropy_array(n) = Q_entropy_volumetric_avg * Volume;
    Q_gen_array(n) = Q_joule_array(n) + Q_entropy_array(n);
    
    % --- Temperature Update Using Finite Differences ---
    for i = 2:Nr-1
        for j = 2:Nz-1
            % Radial conduction term (axisymmetric formulation)
            radial_term = (T(i+1, j) - 2*T(i, j) + T(i-1, j)) / Dr^2 + ...
                          (1/(i*Dr)) * (T(i+1, j) - T(i-1, j)) / (2*Dr);
            % Axial conduction term
            axial_term = (T(i, j+1) - 2*T(i, j) + T(i, j-1)) / Dz^2;
            
            % Local entropy heating term at grid point
            Q_entropy_local = (I * T(i, j) * (delta_S) / 96485) / Volume;
            
            % Total local heat generation per unit mass (W/kg)
            heat_gen = (Q_joule_volumetric + Q_entropy_local) / (rho * c);
            
            % Explicit Euler update for temperature at grid point
            T_new(i, j) = T(i, j) + Dt * (k/(rho*c)) * (radial_term + axial_term) + Dt * heat_gen;
        end
    end
    
    % --- Calculate the RHS of the Heat Equation at the Cell Center ---
    radial_term_center = 4 * (T(2, center_j) - T(1, center_j)) / Dr^2;
    axial_term_center = (T(1, center_j+1) - 2*T(1, center_j) + T(1, center_j-1)) / Dz^2;
    Q_entropy_local_center = (I * T(1, center_j) * (delta_S) / 96485) / Volume;
    Q_gen_center = Q_joule_volumetric + Q_entropy_local_center;
    RHS_term_center = k * (radial_term_center + axial_term_center) + Q_gen_center;
    RHS_term_array(n) = RHS_term_center;
    
    % --- Apply Boundary Conditions ---
    % Axisymmetric boundary at r = 0
    T_new(1, :) = T_new(2, :);
    % Convective (Robin) boundary at r = R
    T_new(end, :) = (k * T(end-1, :) + h * Dr * T_ambient) / (k + h * Dr);
    % Convective boundary conditions at z = 0 and z = L
    T_new(:, 1) = (k * T(:, 2) + h * Dz * T_ambient) / (k + h * Dz);
    T_new(:, end) = (k * T(:, end-1) + h * Dz * T_ambient) / (k + h * Dz);

    
    % Update the temperature field
    T = T_new;
    T_center(n) = T(1, center_j);
    
    % Capture the transition time from CC to CV when SoC reaches the threshold
    if ~found_t1 && SoC >= SOC_threshold
        t1 = current_time;  % Capture the transition time
        found_t1 = true;
        T_t1 = T;         % Store temperature field at t1
    end
    
    % At final time step, store final temperature field
    if n == Nt
        T_final = T;
    end
    
%{

    % Safety check: terminate if temperature exceeds 60°C
    if max(T(:)) > 60
        disp('Temperature exceeded safety limit. Terminating charging.');
        break;
    end

%}

    
    % Visualization update every 50 time steps
    if mod(n, 200) == 0
        [R_mesh, Z_mesh] = meshgrid(linspace(0, R, Nr), linspace(0, L, Nz));
        surf(R_mesh, Z_mesh, T', 'EdgeColor', 'none');
        colorbar;
        xlabel('r (m)'); ylabel('z (m)'); zlabel('Temperature (°C)');
        title(['Time = ', num2str(current_time), ' s, Mode = ', charge_mode]);
        view(2);
        pause(0.1);
    end
end

%% Post-Processing Plots

% Plot Current and Voltage vs. Time
figure;
subplot(2,1,1);
plot(time_array, current_array, 'b', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Current (A)');
title('Current vs. Time');
grid on;
if ~isnan(t1)
    hold on; xline(t1, '--r', ['t_{CV} = ' num2str(t1) ' s']); hold off;
end

subplot(2,1,2);
plot(time_array, voltage_array, 'r', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Voltage (V)');
title('Voltage vs. Time');
grid on;
if ~isnan(t1)
    hold on; xline(t1, '--r', ['t_{CV} = ' num2str(t1) ' s']); hold off;
end

% Plot Heat Generation vs. Time
figure;
plot(time_array, Q_gen_array, 'k', 'LineWidth', 3, 'DisplayName', 'Total Q');
hold on;
plot(time_array, Q_joule_array, 'b--', 'LineWidth', 2, 'DisplayName', 'Joule Heating');
plot(time_array, Q_entropy_array, 'r:', 'LineWidth', 2, 'DisplayName', 'Entropy Heating');
xlabel('Time (s)'); ylabel('Heat Generation (W)');
title('Heat Generation vs. Time');
legend('Location', 'best'); grid on;
if ~isnan(t1)
    xline(t1, '--m', ['t_{CV} = ' num2str(t1) ' s'], 'LineWidth', 1.5);
end
hold off;

% Plot Center Temperature vs. Time
figure;
plot(time_array, T_center, 'm', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Temperature (°C)');
title('Center Temperature vs. Time');
grid on;
if ~isnan(t1)
    hold on; xline(t1, '--r', ['t_{CV} = ' num2str(t1) ' s']); hold off;
end

% Plot RHS Term of the Heat Equation at Cell Center vs. Time
figure;
plot(time_array, RHS_term_array, 'b', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('k(\nabla^2T) + Q_{gen} (W/m^3)');
title('Heat Equation RHS at Cell Center vs. Time');
grid on;
if ~isnan(t1)
    hold on; xline(t1, '--r', ['t_{CV} = ' num2str(t1) ' s'], 'LineWidth', 1.5); hold off;
end

disp(['CC-CV transition occurred at t = ', num2str(t1), ' s']);

%% Temperature Contour Plots (2D)

% Create meshgrid for contour plots
r = linspace(0, R, Nr);
z = linspace(0, L, Nz);
[R_mesh, Z_mesh] = meshgrid(r, z);

figure;

% Temperature Contour at t1 (End of CC Phase)
subplot(1,2,1);
contourf(R_mesh, Z_mesh, T_t1', 20, 'LineColor', 'none');
colorbar;
title(['Temperature at t = ', num2str(t1), ' s (CC to CV Transition)']);
xlabel('Radial Distance (m)');
ylabel('Axial Distance (m)');

% Temperature Contour at t_max (End of Simulation)
subplot(1,2,2);
contourf(R_mesh, Z_mesh, T_final', 20, 'LineColor', 'none');
colorbar;
title(['Temperature at t = ', num2str(total_time), ' s (Final)']);
xlabel('Radial Distance (m)');
ylabel('Axial Distance (m)');

%% Plot SoC vs. Time
figure;
plot(time_array, SoC_array, 'g', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('State of Charge (SoC)');
title('State of Charge vs. Time');
grid on;
if ~isnan(t1)
    hold on; xline(t1, '--r', ['t_{CV} = ' num2str(t1) ' s']); hold off;
end

% Plot Temperature vs. Voltage
figure;
plot(voltage_array,T_center, 'b', 'LineWidth', 3);
xlabel('Battery Voltage (V)');
ylabel('Center Temperature (°C)');
title('Center Temperature vs. Battery Voltage during Charging');
grid on;
hold on;
if ~isnan(t1)
    t1_idx = find(time_array >= t1, 1);
    plot(voltage_array(t1_idx), T_center(t1_idx), 'ko', 'MarkerSize', 5, ...
        'MarkerFaceColor', 'y', 'DisplayName', 'CC-CV Transition');
end
hold off;
legend('Location', 'best');

% Plot Center Temperature vs. SoC
figure;
plot(SoC_array, T_center, 'r', 'LineWidth', 2);
xlabel('State of Charge (SoC)');
ylabel('Center Temperature (°C)');
title('Center Temperature vs. State of Charge during Charging');
grid on;
xlim([min(SoC_array) 1]);
hold on;
if ~isnan(t1)
    t1_idx = find(time_array >= t1, 1);
    plot(SoC_array(t1_idx), T_center(t1_idx), 'ko', 'MarkerSize', 5, ...
        'MarkerFaceColor', 'y', 'DisplayName', 'CC-CV Transition');
end
hold off;
legend('Location', 'best');


%%TEMP MATRIX
%% Export Results to Excel
%{
% Time-series data
results_table = table(time_array', current_array', voltage_array', SoC_array', ...
                      Q_joule_array', Q_entropy_array', Q_gen_array', ...
                      T_center', RHS_term_array', ...
                      'VariableNames', {'Time_s', 'Current_A', 'Voltage_V', ...
                                        'SoC', 'Q_joule_W', 'Q_entropy_W', ...
                                        'Q_total_W', 'T_center_degC', 'RHS_term'});

% Save to Excel file
writetable(results_table, 'battery_thermal_results.xlsx');


if ~isempty(T_t1)
    writematrix(T_t1, 'Temperature_at_t1.csv');
end

if ~isempty(T_final)
    writematrix(T_final, 'Temperature_Final.csv');
end
%}



%temp at each interval (imp)
% Downsample data to 1-second intervals
%{
step_interval = round(1 / Dt);  % 100 if Dt = 0.01
indices = 1:step_interval:length(time_array);

% Create a table for export
output_table = table();
output_table.Time_s = time_array(indices)';
output_table.Current_A = current_array(indices)';
output_table.Voltage_V = voltage_array(indices)';
output_table.SoC = SoC_array(indices)';
output_table.Q_joule_W = Q_joule_array(indices)';
output_table.Q_entropy_W = Q_entropy_array(indices)';
output_table.Q_total_W = Q_gen_array(indices)';
output_table.CenterTemp_C = T_center(indices)';
output_table.RHS_center = RHS_term_array(indices)';

% Export to Excel
filename = 'battery_thermal_results_per_sec.xlsx';
writetable(output_table, filename);
disp(['Results saved to ', filename]);
%}
