clear all
close all
%discharge-single cell
%% Thermal and Geometric Parameters
L = 0.1;         % Cell length (m)
R = 0.01;        % Cell radius (m)
k = 3;           % Thermal conductivity (W/mK)
rho = 2700;      % Density (kg/m^3)
c = 800;         % Specific heat (J/kg.K)
h = 10;          % Convective heat transfer coefficient (W/m^2K)
T_ambient = 25;  % Ambient temperature (°C)
C_rate = 2;

%% Battery Parameters (Same as charging)
E0 = 4.3;          
K = 0.001;         
Q_max = 10;         
R_internal = 0.01; 

I_discharge = Q_max * C_rate;      % Constant discharge current (A)

V_max = 4.2;       % Maximum charging voltage (V)
V_min = 2.7;           % Cutoff voltage for discharge (V)

delta_S = -5;         
F = 96485;             

SoC = 1.0;              % Full charge initially
Q_current = SoC * Q_max;

V_batt = E0 - K * ((Q_max - Q_current) / Q_current) * I_discharge - R_internal * I_discharge;

%% Derived Parameters
Volume = pi * R^2 * L;

Nr = 50;
Nz = 100;
Dr = R / (Nr - 1);
Dz = L / (Nz - 1);
Dt = 0.01;
total_time = 4000;
Nt = total_time / Dt;

T = T_ambient * ones(Nr, Nz); 
T_new = T;

%% Preallocate Arrays
time_array = (0:Nt-1) * Dt;
current_array = zeros(1, Nt);
voltage_array = zeros(1, Nt);
Q_joule_array = zeros(1, Nt);
Q_entropy_array = zeros(1, Nt);
Q_gen_array = zeros(1, Nt);
T_center = zeros(1, Nt);
SoC_array = zeros(1, Nt);

center_j = ceil(Nz/2);
stop_flag = false;

%% Time-Stepping Loop
for n = 1:Nt
    current_time = n * Dt;

    % Update SoC
    SoC = SoC - (I_discharge * Dt / 3600) / Q_max;
    if SoC < 0, SoC = 0; end
    Q_current = SoC * Q_max;

    % Voltage update (Shepherd model)
    if Q_current > 0
        V_batt = E0 - K * ((Q_max - Q_current) / Q_current) * I_discharge - R_internal * I_discharge;
    else
        V_batt = V_min; % Prevent divide-by-zero
    end

    % Cutoff condition
    if V_batt <= V_min
        disp(['Discharge complete at time = ', num2str(current_time), ' s']);
        stop_idx = n;
        stop_flag = true;
    end

    current_array(n) = I_discharge;
    voltage_array(n) = V_batt;
    SoC_array(n) = SoC;

    % Heat generation
    T_avg = mean(T(:));
    Q_joule_volumetric = (I_discharge^2 * R_internal) / Volume;
    Q_entropy_volumetric_avg = (I_discharge * T_avg * delta_S / F) / Volume;

    Q_joule_array(n) = Q_joule_volumetric * Volume;
    Q_entropy_array(n) = Q_entropy_volumetric_avg * Volume;
    Q_gen_array(n) = Q_joule_array(n) + Q_entropy_array(n);

    % --- Temperature Update ---
    for i = 2:Nr-1
        for j = 2:Nz-1
            radial_term = (T(i+1, j) - 2*T(i, j) + T(i-1, j)) / Dr^2 + ...
                          (1/(i*Dr)) * (T(i+1, j) - T(i-1, j)) / (2*Dr);
            axial_term = (T(i, j+1) - 2*T(i, j) + T(i, j-1)) / Dz^2;
            Q_entropy_local = (I_discharge * T(i, j) * delta_S / F) / Volume;
            heat_gen = (Q_joule_volumetric + Q_entropy_local) / (rho * c);

            T_new(i,j) = T(i,j) + Dt * (k/(rho*c)) * (radial_term + axial_term) + Dt * heat_gen;
        end
    end

    % --- Boundary Conditions ---
    T_new(1, :) = T_new(2, :);
    T_new(end, :) = (k * T(end-1, :) + h * Dr * T_ambient) / (k + h * Dr);
    T_new(:, 1) = (k * T(:, 2) + h * Dz * T_ambient) / (k + h * Dz);
    T_new(:, end) = (k * T(:, end-1) + h * Dz * T_ambient) / (k + h * Dz);

    T = T_new;
    T_center(n) = T(1, center_j);

    % Visualization every 200 steps
    if mod(n, 200) == 0
        [R_mesh, Z_mesh] = meshgrid(linspace(0, R, Nr), linspace(0, L, Nz));
        surf(R_mesh, Z_mesh, T', 'EdgeColor', 'none');
        xlabel('r (m)'); ylabel('z (m)'); zlabel('Temperature (°C)');
        title(['Time = ', num2str(current_time), ' s']);
        view(2); colorbar;
        pause(0.1);
    end

    if stop_flag
        break
    end
end

%% Post-Processing
figure;
subplot(2,1,1);
plot(time_array(1:n), current_array(1:n), 'b', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Current (A)');
title('Constant Current Discharge'); grid on;

subplot(2,1,2);
plot(time_array(1:n), voltage_array(1:n), 'r', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Voltage (V)');
title('Battery Voltage During Discharge'); grid on;

figure;
plot(time_array(1:n), Q_gen_array(1:n), 'k', 'LineWidth', 3);
hold on;
plot(time_array(1:n), Q_joule_array(1:n), 'b--', 'LineWidth', 2);
plot(time_array(1:n), Q_entropy_array(1:n), 'r:', 'LineWidth', 2);
legend('Total Q', 'Joule', 'Entropy'); xlabel('Time (s)'); ylabel('Heat (W)');
title('Heat Generation During Discharge'); grid on;

figure;
plot(time_array(1:n), T_center(1:n), 'm', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Center Temp (°C)');
title('Center Temperature vs Time During Discharge'); grid on;

figure;
plot(voltage_array(1:n), T_center(1:n), 'b', 'LineWidth', 3);
xlabel('Battery Voltage (V)');
ylabel('Center Temperature (°C)');
title('Center Temperature vs Battery Voltage During Discharge');
grid on;

% Add this line to reverse the x-axis direction
set(gca, 'XDir', 'reverse');

%% Plot SoC vs. Time
figure;
plot(time_array(1:n), SoC_array(1:n), 'g', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('State of Charge (SoC)');
title('State of Charge vs. Time');
grid on;

% Plot Center Temperature vs. SoC
figure;
plot(SoC_array(1:n), T_center(1:n), 'r', 'LineWidth', 2);
xlabel('State of Charge (SoC)');
ylabel('Center Temperature (°C)');
title('Center Temperature vs SoC During Discharge');
grid on;
% Add this line to reverse the x-axis direction (since SoC decreases during discharge)
set(gca, 'XDir', 'reverse');
