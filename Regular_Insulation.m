clc; clearvars;

TimeLoc = 'AUSTIN 08.17.2023';
parts = split(TimeLoc, ' ');
city = parts{1};
dateStr = parts{2};

TempSpreadsheetNew = getCityWeatherData(city, dateStr, 'TemperatureDataCity.xlsx');

% Given parameters
A = 10; % Area in m^2

% Limestone (left side ) Properties
k1 = 1.3; % Thermal conductivity (W/m·K)
cp_k1 = 0.91; % Specific heat (kJ/kg·K) 
roh_k1 = 1550; % Density (kg/m^3)

% Drywall (right side) Properties
k2 = 0.21; % Thermal conductivity (W/m·K)
cp_k2 = 0.95; % Specific heat (kJ/kg·K)
roh_k2 = 700; % Density (kg/m^3)

% Insulation properties
kinsulation = 0.04; % Thermal conductivity (W/m·K)
cp_insulation = 0.84; % Specific heat (kJ/kg·K)
roh_insulation = 12; % Density (kg/m^3)

% Thicknesses (m)
L1 = 0.1; % Limestone
Linsulation = 0.012; % thickness insulation
L2 = 0.0127; % Drywall thickness

Mk1 = roh_k1 * L1; % kG
Mk2 = roh_k2 * L2; % kG
Mk3 = roh_insulation * Linsulation;

% Temperatures
T_inside = 20; % °C 
T_outside = TempSpreadsheetNew{1, 6}; % °C
T_mid = (T_inside + T_outside) / 2; % Inital mid temp

% Time Parameters
dt = 1; % Time step in seconds
t_max = 60*60*24; % Total time (seconds)
time_steps = t_max / dt; % Number of time steps

% Define resistance function
getResistance = @(L, k, A) L / (k * A);
mCAT = @(Q, m, c) Q / (m * c);

% flags
Temps = zeros(1, time_steps);
MeltFraction = 0;
tol = 1e-2;
counter = 1;
vals = zeros(60*60*24, 1);

RLimestone = getResistance(L1, k1, A); % returns in units of seconds * kelvin  / Joules (constant)
RDrywall = getResistance(L2, k2, A); % returns in units of seconds * kelvin  / Joules (constant)
RInsulation = getResistance(Linsulation, kinsulation, A); 

% Calculating Heat flux
q_net = zeros(1, time_steps);

% Build a time array for the outdoor temp data
nWeatherPts = height(TempSpreadsheetNew);
time_weather = linspace(0, t_max, nWeatherPts); % Seconds
temp_weather = TempSpreadsheetNew{:, 6};

% Interpolate T_outside for every second
T_outside_all = interp1(time_weather, temp_weather, (1:time_steps));

for t = 1:time_steps

        T_outside = T_outside_all(t);
    
        Q_in = (T_outside - T_mid) / ((RInsulation/2) + RLimestone); % J/s
        Q_out = (T_mid - T_inside) / (RDrywall+(RInsulation/2));           % J/s
        Q_net = Q_in - Q_out;  % net energy
        T_mid = T_mid + ( (Q_net / 1000) / ((Mk1 * cp_k1) + (Mk2 * cp_k2) + (Mk3 * cp_insulation)) );

        Temps(t) = T_mid; % Tracks Temps 

        q_net(t) = (Q_in / A) - (Q_out / A);  % or just store this separately if you prefer


end

avgPCMTemp = mean(Temps);
maxRateChange = max(abs(diff(Temps))) * 3600; % °C/hr

% Display variables
dateFormatted = strrep(dateStr, '.', '/');

figure;
tiledlayout(2,2);
sgtitle(['Simulation for ' upper(city) ' on ' dateFormatted ' (Regular Insulation)']);

% Plot 2
nexttile;
plot((1:time_steps)/3600, Temps, 'b', 'LineWidth', 2);
title('Midpoint Temperature');

% Plot 3
nexttile;
plot((1:time_steps)/3600, T_outside_all, 'b', 'LineWidth', 2);
title('Outside Temperature');

q_net_smooth = movmean(q_net, 60); % 60-second smoothing window

% Plot 4 - Heat Flux
nexttile;
plot((1:time_steps)/3600, q_net_smooth, 'b', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Heat Flux (W/m²)');
title('Heat Flux Through Wall Layers');

% Summary Stats'
nexttile;
axis off;
y = 0.8;
dy = 0.17;

text(0, y - 2*dy, sprintf('Avg Insulation Temp: %.2f °C', avgPCMTemp));
text(0, y - 4*dy, sprintf('Max Temp Rate Change: %.2f °C/hr', maxRateChange));
text(0, y - 3*dy, sprintf('Max Insulation Temp: %.2f °C', max(Temps)));
title('Summary Stats');