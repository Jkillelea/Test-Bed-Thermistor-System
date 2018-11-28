% Find the optimimum second resistor for a 10K thermistor
clear all; clc;

Vin = 3.3;
Tmax =  60;
Tmin = -30;

adc_bits = 12;
adc_bins = 2^adc_bits - 1;
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
voltage_range = min_max_temp_reading(:, 1) - min_max_temp_reading(:, 2);
[maxval, idx] = max(voltage_range);
ideal_resistance = resistance_range(idx);

fprintf('Ideal resistance %d kOhms\n', ideal_resistance/1000);

% Try and plot only if Java UI libs are present
if usejava('awt') | usejava('swing')
    figure; hold on; grid on;
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
end

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

figure; hold on; grid on;
title('Voltage sensitivty to temperature');
plot(temps(2:end), dVdC(2:end))

figure; hold on; grid on;
title('ADC bins per deg C vs temperature');
plot(temps(2:end), bins_per_degC(2:end), 'linewidth', 2)

figure; hold on; grid on;
title('Deg C per bin vs temperature');
plot(temps(2:end), degC_per_bin(2:end), 'linewidth', 2)
