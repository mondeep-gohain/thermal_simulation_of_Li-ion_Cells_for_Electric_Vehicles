clear all
close all

%% Thermal and Geometric Parameters
L = 0.065;         % Length of the cylindrical cell (m)
D = 0.018;       % Diameter of the cell (m)
R = D/2;        % Radius of the cylindrical cell (m)
k = 3;           % Thermal conductivity (W/mK)
rho = 2871.74;      % density (mass/volume)
c = 859;         % Specific heat (J/kg.K)
h = 10;          % Convective heat transfer coefficient (W/m^2K)
T_ambient = 27;  % Ambient temperature (°C)

%% Battery Parameters
E0 = 3.6;          % Constant voltage term (V)
K = 0.001;         % Polarization constant (V/Ah)
Q_cell = 2.9;       % Capacity of one cell (Ah)
R_internal = 0.02; % Internal resistance (Ohm)
V_max = 4.2;       % Maximum voltage per cell (V)
delta_S = -5;    % Entropy change (J/mol.K) ~ tuned
F = 96485;         % Faraday’s constant (C/mol)
tau_CV = 50;       % Time constant for CV phase (s)
SOC_threshold = 0.85;
C_rate = 0.5;

%% Pack configuration
Ns = 4;  % Series cells
Np = 2;  % Parallel strings
Ncells = Ns * Np;

Q_pack = Np * Q_cell;      % Pack capacity (Ah)
I_total_CC = Q_pack * C_rate;   % 1C current for pack = 40 A
I_min = Q_pack * (1/10);      % Minimum CV current (A)

%% Initial conditions
SoC = 0.2;   % all cells start equal

%% Derived parameters
Volume = pi * R^2 * L; % Volume of one cell (m^3)

% Time stepping
Dt = 0.01;         
total_time = 10000; 
Nt = total_time / Dt;

%% Initialize arrays
time_array = (0:Nt-1) * Dt;
current_array = zeros(1, Nt);   % Pack current
voltage_array = zeros(1, Nt);   % Pack voltage
SoC_array = zeros(1, Nt);       % SOC (per cell, all same)
Q_gen_vol_array = zeros(1, Nt); % Heat per cell

found_t1 = false;
t1 = NaN;

%% Simulation
I_total = I_total_CC; % Pack current in CC phase

for n = 1:Nt
    current_time = n * Dt;

        % Stop simulation if fully charged
    if SoC >= 1.0 %   && I_total <= I_min
        Nt = n; % Truncate arrays later
        break;
    end

    if SoC < SOC_threshold
        % --- CC Phase ---
        I_pack = I_total;          % pack current (40 A at 1C)
        I_string = I_pack / Np;    % each string current
        I_cell = I_string;         % series cells: same current (10 A)

        SoC = SoC + (I_cell * Dt / 3600) / Q_cell;
        if SoC > 1, SoC = 1; end

        Q_current = SoC * Q_cell;
        V_cell = E0 - K * ((Q_cell - Q_current) / Q_current) * I_cell - R_internal * I_cell;
        if V_cell > V_max, V_cell = V_max; end
        V_pack = Ns * V_cell;

    else
        % --- CV Phase ---
        V_cell = V_max;
        V_pack = Ns * V_cell;
        I_total = max(I_total * exp(-Dt / tau_CV), I_min); % pack current decays
        I_pack = I_total;
        I_string = I_pack / Np;
        I_cell = I_string;

        SoC = SoC + (I_cell * Dt / 3600) / Q_cell;
        if SoC > 1, SoC = 1; end
    end

    current_array(n) = I_pack;
    voltage_array(n) = V_pack;
    SoC_array(n) = SoC;

    % Heat generation per cell
    Q_joule_vol = (I_cell^2 * R_internal) / Volume;
    Q_entropy_vol = (I_cell * (T_ambient+273.15) * delta_S / F) / Volume;
    Q_gen_vol = Q_joule_vol + Q_entropy_vol;
    Q_gen_vol_array(n) = Q_gen_vol;

    if ~found_t1 && SoC >= SOC_threshold
        t1 = current_time;
        found_t1 = true;
    end
end

%% Save heat generation (per cell)
step_interval = round(1 / Dt);
indices = 1:step_interval:Nt;
time_array_ansys = time_array(indices);
Q_gen_vol_ansys = Q_gen_vol_array(indices);

data = [time_array_ansys', Q_gen_vol_ansys'];
writematrix(data, 'heat_generation_NCA.csv');
disp('Heat generation per cell saved to heat_generation_NCA.csv');

%% Save SOC vs time with 1-second intervals (up to same end time as heat table)
time_step_ansys = 1;  % match ANSYS output
end_time = time_array_ansys(end);  % same end time as heat generation data
time_resampled = (0:time_step_ansys:end_time)';  % only till that time

soc_resampled = interp1(time_array(1:Nt), SoC_array(1:Nt), time_resampled, 'linear');
soc_data = [time_resampled, soc_resampled];

writematrix(soc_data, 'soc_vs_time_NCA_charging.csv');
disp('SOC vs time saved (trimmed to same end time as heat generation table).');

