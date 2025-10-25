% List of files
files = {
    'Temp_vs_SoC_NMC_0p5C_3S.csv'
    'Temp_vs_SoC_NMC_1C_3S.csv'
    'Temp_vs_SoC_NMC_1p5C_3S.csv'
    'Temp_vs_SoC_NMC_2p5C_3S.csv'
    'Temp_vs_SoC_NMC_5C_3S.csv'
};

C_rates = {'0.5C','1C','1.5C','2.5C','5C'};   % labels for each file

% Read first file to get SoC reference
data0 = readmatrix(files{1});  
SoC_ref = data0(:,1);          
Combined = table(SoC_ref,'VariableNames',{'SoC'});  

% Loop over all files
for i = 1:numel(files)
    data = readmatrix(files{i});
    SoC = data(:,1);
    T   = data(:,2);

    % Remove rows with NaN or Inf
    valid = isfinite(SoC) & isfinite(T);
    SoC = SoC(valid);
    T   = T(valid);

    % Ensure unique & sorted SoC
    [SoC, idx] = unique(SoC,'stable');
    T = T(idx);

    % Interpolate temperature onto reference SoC grid
    T_interp = interp1(SoC, T, SoC_ref, 'linear', 'extrap');
    
    % Add as new column
    Combined.(C_rates{i}) = T_interp;
end

% Save combined CSV
writetable(Combined,'Combined_Temp_vs_SoC_NMC_charge_3S.csv');
disp('✅ Combined table saved as Combined_Temp_vs_SoC.csv');
