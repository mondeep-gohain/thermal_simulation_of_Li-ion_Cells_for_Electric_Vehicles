clear all
close all

%% Thermal and Geometric Parameters
L = 0.065;         % Length of the cylindrical cell (m)
D = 0.018;       % Diameter of the cell (m)
R = D/2;        % Radius of the cylindrical cell (m)
k = 3;           % Thermal conductivity (W/mK)
rho = 2962.42;      % density (mass/volume)
c = 997;         % Specific heat (J/kg.K)
h = 10;          % Convective heat transfer coefficient (W/m^2K)
T_ambient = 27;  % Ambient temperature (°C)

%% Battery Parameters
E0 = 3.635;          % Constant voltage term (V)
K = 0.001;         % Polarization constant (V/Ah)
Q_cell = 3.5;       % Capacity of one cell (Ah)
R_internal = 0.02; % Internal resistance (Ohm)
V_min = 2.7;       % Cutoff voltage per cell (V)
delta_S = -5;    % Entropy change (J/mol.K) ~ tuned
F = 96485;         % Faraday’s constant (C/mol)
C_rate = 2.5;

%% Pack configuration
Ns = 4;  % 4 cells in series
Np = 1;  % single string
Ncells = Ns * Np;

Q_pack = Np * Q_cell;  % Pack capacity (Ah)
I_discharge = Q_pack * C_rate;      % 2C discharge current per cell (A)

%% Initial conditions
SoC = 1.0;   % start fully charged

%% Derived parameters
Volume = pi * R^2 * L; % Volume of one cell (m^3)

% Time stepping
Dt = 0.1;         
total_time = 8000; 
Nt = total_time / Dt;

%% Initialize arrays
time_array = (0:Nt-1) * Dt;
current_array = zeros(1, Nt);   % Pack current
voltage_array = zeros(1, Nt);   % Pack voltage
SoC_array = zeros(1, Nt);       % SOC
Q_gen_vol_array = zeros(1, Nt); % Heat per cell

%% Discharge Simulation
for n = 1:Nt
    current_time = n * Dt;

    I_cell = I_discharge;  % same current through each cell
    SoC = SoC - (I_cell * Dt / 3600) / Q_cell;
    if SoC < 0, SoC = 0; end

    Q_current = SoC * Q_cell;

    % Cell voltage (Shepherd discharge model)
    if Q_current > 0
        V_cell = E0 - K * (Q_cell / (Q_current)) * I_cell - R_internal * I_cell;
    else
        V_cell = V_min; % prevent numerical blow-up
    end

    % Pack voltage
    V_pack = Ns * V_cell;

    % Save data
    current_array(n) = I_cell;
    voltage_array(n) = V_pack;
    SoC_array(n) = SoC;

    % Heat generation per cell
    Q_joule_vol = (I_cell^2 * R_internal) / Volume;
    Q_entropy_vol = (I_cell * (T_ambient+273.15) * delta_S / F) / Volume;
    Q_gen_vol = Q_joule_vol + Q_entropy_vol;
    Q_gen_vol_array(n) = Q_gen_vol;

    % Stop if cell hits cutoff voltage
    if V_cell <= V_min
        break;
    end
end

%% Save Heat Generation Table for ANSYS
step_interval = round(1 / Dt);
valid_indices = 1:step_interval:n; % only keep simulated steps
time_array_ansys = time_array(valid_indices);
Q_gen_vol_ansys = Q_gen_vol_array(valid_indices);

data = [time_array_ansys', Q_gen_vol_ansys'];
writematrix(data, 'heat_generation_NMC_discharge_.csv');
disp('Heat generation per cell during discharge saved to heat_generation_nmc_discharge_.csv');

%% Save SOC vs time with 1-second intervals
time_step_ansys = 1;  % match ANSYS output
time_resampled = (0:time_step_ansys:max(time_array))';

soc_resampled = interp1(time_array, SoC_array, time_resampled, 'linear');

soc_data = [time_resampled, soc_resampled];
writematrix(soc_data, 'soc_vs_time_NMC_discharging.csv');
disp('SOC vs time');
