% Find the optimimum second resistor for a 10K thermistor
clear all; clc; close all;

Vin           = 3.3;
Tmax          =  60;
Tmin          = -30;
adc_bits      = 12;
adc_bins      = 2^adc_bits - 1;
volts_per_bin = Vin / adc_bins;

R25 = 10e3; % Ohms
B25 = 3435; % Kelvin
T25 = 273.15 + 25; % 25 degC in Kelvin

% Load data between Tmin and Tmax
data = load('./thermistor_spec');
idx  = Tmin <= data(:, 1) & data(:, 1) <= Tmax;
data = data(Tmin <= data(:, 1) & data(:, 1) <= Tmax, :);
resistance_range = linspace(1e3, 1e5, 1e4);

min_max_temp_reading = zeros(length(resistance_range), 2);

for i = 1:length(resistance_range)
    R1 = resistance_range(i);
    
    % -30 degC
    Rth = data(1, 2);
    output_voltage = Vin*(Rth/(Rth+R1));
    min_max_temp_reading(i, 1) = output_voltage;

    % 60 degC
    Rth = data(end, 2);
    output_voltage = Vin*(Rth/(Rth+R1));
    min_max_temp_reading(i, 2) = output_voltage;
end

% Find the resistor R1 that gives the widest voltage range across the relevant temperatures
voltage_range    = min_max_temp_reading(:, 1) - min_max_temp_reading(:, 2);
[maxval, idx]    = max(voltage_range);
ideal_resistance = resistance_range(idx);

fprintf('Ideal resistance %.2f kOhms\n', ideal_resistance/1000);

f = figure; hold on; grid on;
title('Voltage Range vs Resistor Value');
xlabel('log_{10}(R1) Resistance (log_{10}(\Omega))');
ylabel('Output Voltage');
semilogx(resistance_range, min_max_temp_reading(:, 1), 'b', ... 
                        'displayname', 'Voltage at -30 degC')
semilogx(resistance_range, min_max_temp_reading(:, 2), 'r', ... 
                        'displayname', 'Voltage at 60 degC')
semilogx(resistance_range, voltage_range, 'g', ... 
                        'linewidth', 2, ... 
                        'displayname', 'Voltage Range')
semilogx(ideal_resistance, maxval, 'ro')
legend('show')
print(f, 'images/voltage_range_vs_resistor_value', '-dpng')

% Calculate resolutions
data          = load('./thermistor_spec');
R1            = 18e3;                     % Chosen resistor spec
Rth           = data(:, 2);               % Thermistor R
temps         = data(:, 1);               % Temperatures at which R is measured
Vout          = Vin .* Rth ./ (Rth + R1); % Vout at each temperature
bins_per_volt = 1/volts_per_bin;          % Voltage resolution of the ADC

% calculate senitivity of voltage to temperature
dVdC = zeros(length(Vout), 1);
for i = 2:length(Vout)
    dV = abs(Vout(i) - Vout(i-1));
    dC = abs(temps(i) - temps(i-1));
    dVdC(i) = dV/dC;
end

% Temperature resolution
bins_per_degC = bins_per_volt .* dVdC;
degC_per_bin = 1 ./ bins_per_degC;

f = figure; hold on; grid on;
title('Voltage sensitivty to temperature');
plot(temps(2:end), dVdC(2:end), 'linewidth', 2)
print(f, 'images/voltage_sensitivity_to_temperature', '-dpng')

f = figure; hold on; grid on;
title('ADC bins per deg C vs temperature');
title('ADC bins');
title('Degrees C');
plot(temps(2:end), bins_per_degC(2:end), 'linewidth', 2)
print(f, 'images/adc_bins_per_deg_C_vs_temperature', '-dpng')

f = figure; hold on; grid on;
title('Deg C per bin vs temperature');
plot(temps(2:end), degC_per_bin(2:end), 'linewidth', 2)
print(f, 'images/deg_c_per_bin_vs_temperature', '-dpng')

% Error analysis
%% N = number of samples, t = coeff for 95% condfidence interval
N = 1; t = 12.71;
% N = 3; t = 3.128;
% N = 10; t = 2.228;
% N = 100; t = 1.984;

Rth_bias   = Rth .* (data(:, 4)./100); % Deviation from correct (percent err times ohms)
R1_error   = R1 * 0.01; % 1% error
Vin_error  = 0.001;     % 1mV? TODO
Vout_error = 0.004 .* Vout; % (-0.6% to +0.1%) + 1 mV noise
Vout_noise = 0.001;     % 1 mV
B25_error  = 0.01 * B25; % 1%, Kelvin
R25_error  = 0.01 * R25; % 1%, Ohms

dRthdVin  = -R1 ./ (Vout .* (Vin./Vout - 1).^2);
dRthdVout = R1*Vin ./ ((Vin./Vout - 1).^2 .* Vout.^2);
dRthdR1   = (Vin./Vout - 1).^-1;

% Rth uncertainty
Rth_bias_err = sqrt((dRthdVin  .* Vin_error).^2 ...
                  + (dRthdVout .* Vout_error).^2 ...
                  + (dRthdR1   .* R1_error).^2) ...
                  + abs(Rth_bias);

Rth_rand_err = abs(dRthdVout .* Vout_noise);

Rth_Unc = t*Rth_rand_err/sqrt(N) + Rth_bias_err;

f = figure; hold on; grid on;
title(sprintf('Rth uncertainty vs Rth (normalized), N=%d', N))
xlabel('Thermistor Resistance (Ohms)');
ylabel('Thermistor Uncertainty (Percent)');
plot(Rth, 100*Rth_Unc./Rth, 'linewidth', 2)
print(f, sprintf('images/rth_uncertainty_vs_rth_normalized_N%d', N), '-dpng')


% Temperature error
% Emperical fit
T = @(Rth) (1/T25 + log(Rth/R25)/B25).^-1 - 273.15;

% Sensitivities
dTdB25 = log(Rth/R25) ./ (B25^2 * (1/T25 + log(Rth/R25)/B25).^2);
dTdR25 = 1            ./ (B25*R25*(1/T25 + log(Rth/R25)/B25).^2);
dTdRth = -1           ./ (B25*R25*(1/T25 + log(Rth/R25)/B25).^2);

% Temperature Uncertainty
T_bias_err = sqrt((dTdB25 .* B25_error).^2 ...
                + (dTdR25 .* R25_error).^2 ...
                + (dTdRth .* Rth_bias_err  ).^2);

T_rand_err = abs(dTdRth .* Rth_rand_err);

T_Unc = t*T_rand_err/sqrt(N) + T_bias_err;

f = figure; hold on; grid on;
title(sprintf('Temperature uncertainty vs Temperature, N=%d', N))
xlabel('Temperature (deg C)');
ylabel('Temperature Uncertainty (deg C)');
plot(temps, T_Unc, 'linewidth', 2)
plot([Tmin Tmin], [min(T_Unc) max(T_Unc)], 'r');
plot([Tmax Tmax], [min(T_Unc) max(T_Unc)], 'r');
print(f, sprintf('images/temp_uncertainty_vs_temp_N%d', N), '-dpng')


% Sanity checks
assert(N > 0);
assert(all(Rth_bias_err > 0), 'Rth_bias < 0');
assert(all(Rth_rand_err > 0), 'Rth_rand < 0');
assert(all(Rth_Unc > 0), 'Rth_Unc < 0');
assert(all(T_bias_err > 0), 'T_bias < 0');
assert(all(T_rand_err > 0), 'T_rand < 0');
assert(all(T_Unc > 0), 'T_Unc < 0');
