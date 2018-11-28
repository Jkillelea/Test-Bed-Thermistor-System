% Find the optimimum second resistor for a 10K thermistor
clear all; clc; close all;

Vin           = 3.3;
Tmax          =  60;
Tmin          = -30;
adc_bits      = 12;
adc_bins      = 2^adc_bits - 1;
volts_per_bin = Vin / adc_bins;

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

% Try and plot only if Java UI libs are present
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

% Calculate resolution
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
plot(temps(2:end), bins_per_degC(2:end), 'linewidth', 2)
print(f, 'images/adc_bins_per_deg_C_vs_temperature', '-dpng')

f = figure; hold on; grid on;
title('Deg C per bin vs temperature');
plot(temps(2:end), degC_per_bin(2:end), 'linewidth', 2)
print(f, 'images/deg_c_per_bin_vs_temperature', '-dpng')

% Error analysis
Rth_error = Rth .* (data(:, 4)./100);
R1_error  = R1 * 0.01; % 1% error
Vin_error = 0.001;     % 1mV? TODO

dVoutdVin = Rth ./ (Rth + R1);
dVoutdRth = Vin ./ (Rth + R1) - Vin .* Rth ./ ((Rth + R1).^2);
dVoutdR1  = - Rth .* Vin ./ ((Rth + R1).^2);

f = figure; hold on; grid on;
subplot(1, 2, 1);
title('Rth percent error over temperature');
xlabel('temperature');
ylabel('percent error');
plot(temps, data(:, 4));
subplot(1, 2, 2);
plot(temps, Rth_error);
title('Rth percent error over temperature');
xlabel('temperature');
ylabel('Rth error (ohms)');

f = figure; hold on; grid on;
subplot(2, 2, 1);
plot(temps, dVoutdRth)
subplot(2, 2, 2);
plot(temps, Rth_error)
subplot(2, 2, 3);
plot(temps, dVoutdRth .* Rth_error)


f = figure; hold on; grid on;
title('Rth over temperature')
xlabel('Temperature');
ylabel('Rth');
errorbar(temps, Rth, Rth_error, 'linewidth', 2);
print(f, 'images/rth_over_temperature', '-dpng');

% RMS of sensitivity to terms times uncertainty in each term
Vout_uncertainty = sqrt((dVoutdVin .* Vin_error).^2  ...
                      + (dVoutdR1  .* R1_error ).^2  ...
                      + (dVoutdRth .* Rth_error).^2);

f = figure; hold on; grid on;
title('Vout Uncertainty over Temperature');
xlabel('Temperature');
ylabel('Vout Uncertainty');
plot(temps, Vout_uncertainty, 'linewidth', 2)
print(f, 'images/vout_uncertainty', '-dpng');

p = polyfit(temps, Vout, 4);
g = @(x) p(1) * x.^4 + p(2) * x.^3 + p(3) * x.^2 + p(4) * x + p(5);

f = figure; hold on; grid on;
title('Vout over Temperature curve fit');
xlabel('Temperature');
ylabel('Vout');
plot(temps, Vout, 'linewidth', 2)
plot(temps, g(temps), 'linewidth', 2)
print(f, 'images/vout_vs_temp_curve', '-dpng');

