clc; clearvars;

% The user has to put in the correct file location for PCM_Running_List.xlsx and
% TemperatureDataCity.xlsx.

if exist('C:\Users\pasca\OneDrive\School Work\Spring25\Research\Matlab Programs\PCM_Running_List.xlsx', 'file') ~= 2
    error('Error: file could not be found. Check that the file pathway exists.');
end

if exist('C:\Users\pasca\OneDrive\School Work\Spring25\Research\Matlab Programs\TemperatureDataCity.xlsx', 'file') ~= 2
    error('Error: file could not be found. Check that the file pathway exists.');
end

data = getData('PCM_Running_List.xlsx', 'PCMdata');
PCMname = input('Name of the PCM: ', 's');

cln = find(strcmpi(data(1,:), PCMname), 1);

if isempty(cln)
    error('Error: PCM not found. Please check input or data and retry.');
end

TimeLoc = input('City and Date (in "CITYNAME XX.XX.XXXX" format): ', 's');
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
L1 = 0.1; % length (m)

% Drywall (right side) Properties
k2 = 0.21; % Thermal conductivity (W/m·K)
cp_k2 = 0.95; % Specific heat (kJ/kg·K)
roh_k2 = 700; % Density (kg/m^3)
L2 = 0.0127; % length (m)

% PCM Properties (imported from PCM_Running_List.xlsx location)
Cp_solid = data{2, cln}; % kJ/kg·K 
T_melt = data{3, cln}; % °C
H_fusion = data{4, cln}; % kJ/kg 
Cp_liquid = data{5, cln} ; % kJ/kg·K
k_solid = data{6, cln}; % W/m·K
k_liquid = data{7, cln}; % W/m·K
rho_solid = data{8, cln}; % kg/m^3
rho_liquid = data{9, cln}; % kg/m^3

Dxs = 0.012; % thickness of solid PCM (12 mm, given in source)
Dxl = 0; %thickness of liquid PCM
Dxt = Dxl + Dxs; %total length of PCM

% mass calculations
Msolid = Dxs * A * rho_solid; % mass in KG
Mliquid = Dxl * A * rho_liquid; % mass in KG
Mtotal = Msolid + Mliquid; % mass in KG

Mk1 = roh_k1 * L1; % kG
Mk2 = roh_k2 * L2; % kG

% Temperatures
T_inside = 20; % °C 
T_outside = TempSpreadsheetNew{1, 6}; % °C
T_pcm = 24; % Initial PCM temperature (just below melting value)

% Time Parameters
dt = 1; % Time step in seconds
t_max = 60*60*24; % Total time (seconds)
time_steps = t_max / dt; % Number of time steps

% Define functions
getResistance = @(L, k, A) L / (k * A); % basic resistance function
mCAT = @(Q, m, c) Q / (m * c); % Basic heat absorbtion function

% Initalization for variables used inside model loop
MassPercent = zeros(1, time_steps);
q_net = zeros(1, time_steps);
PCMTemps = zeros(1, time_steps);
OutsideTempMatrix = zeros(60*60*24, 1);
Qleftgraph = zeros(1, time_steps);
Qrightgraph = zeros(1, time_steps);

MeltFraction = 0;
NetQstore = 0;
counter = 1;
verysmall = 0.001;

% Resistance calculations for other materials
RLimestone = getResistance(L1, k1, A); % returns in units of seconds * kelvin  / Joules (constant)
RDrywall = getResistance(L2, k2, A); % returns in units of seconds * kelvin  / Joules (constant)

for t = 1:time_steps

   if counter + 1 <= height(TempSpreadsheetNew)
       if( mod(t, 300) == 0)
            counter = counter + 1;
            T_outside = TempSpreadsheetNew{counter, 6};
       end
       
   end

   RpcmS = getResistance(Dxs, k_solid, A); % returns resistance of solid pcm section in units of seconds * kelvin  / Joules
   RpcmL = getResistance(Dxl, k_liquid, A); % returns resistance of liquid pcm section in units of seconds * kelvin  / Joules
   
   if T_pcm <= T_melt && abs(Mliquid) < verysmall
       % PCM solid and acts as a solid with no phase change while
       % temperature increases
       
        Qleft = (T_outside - T_pcm) / (RpcmL + RLimestone); % J/s
        Qright = (T_pcm - T_inside) / (RpcmS + RDrywall);   % J/s
        Qstore = Qleft - Qright;  % net energy stored in PCM

       T_pcm = T_pcm + ( (Qstore / 1000) / ((Mk1 * cp_k1) + (Mk2 * cp_k2) + (Mtotal * Cp_solid)) );
       
   elseif T_pcm >= T_melt && abs(Msolid) < verysmall
       % PCM liquid and acts as a liquid with no phase change while
       % tempature increases
       
        Qleft = (T_outside - T_pcm) / (RpcmL + RLimestone); 
        Qright = (T_pcm - T_inside) / (RpcmS + RDrywall);           
        Qstore = Qleft - Qright;

        T_pcm = T_pcm + ( (Qstore / 1000) / ((Mk1 * cp_k1) + (Mk2 * cp_k2) + (Mtotal * Cp_liquid)) );


   elseif abs(T_pcm - T_melt) < verysmall
    % PCM undergoing phase change - constant T
       
    Qleft = ((T_outside - T_pcm)) / (RpcmL + RLimestone);
    Qright = ((T_pcm - T_inside)) / (RpcmS + RDrywall);
    Qstore = Qleft - Qright;

    if Qstore >= 0
        % Melting
        Mmelt = Qstore / (H_fusion*1000); % kg
        Mliquid = Mliquid + Mmelt;
        Msolid = Mtotal - Mliquid;
        NetQstore = NetQstore + Qstore;
    else
        % Freezing
        Mfreeze = abs(Qstore) / (H_fusion * 1000); % kg
        Mliquid = Mliquid - Mfreeze;
        Msolid = Mtotal - Mliquid;
    end

    % Clamp values to physical bounds
    Mliquid = max(0, min(Mliquid, Mtotal));
    Msolid = max(0, min(Msolid, Mtotal));

    % Update geometry
    Dxl = (Mliquid / Mtotal) * Dxt;
    Dxs = Dxt - Dxl;
    
    % Kicks out of this loop if PCM is not in phase change by changing the
    % tempature
    if Msolid == 0
        
        % Fully melted → resume as liquid
        Qstore = Qleft - Qright;
        T_pcm = T_pcm + ( (Qstore / 1000) / ((Mk1 * cp_k1) + (Mk2 * cp_k2) + (Mtotal * Cp_liquid)) );
        
    elseif Mliquid == 0
        
        % Fully frozen → resume as solid
        Qstore = Qleft - Qright;
        T_pcm = T_pcm + ( (Qstore / 1000) / ((Mk1 * cp_k1) + (Mk2 * cp_k2) + (Mtotal * Cp_solid)) );
        
    end

   else 
        error('Error: PCM type not accounted for. Check logic.');
   
   end

    % Tracking Arrays
    PCMTemps(t) = T_pcm; % Tracks PCM Temperatures
    MassPercent(t) = Mliquid / Mtotal; % Tracks Mass Fraction
    OutsideTempMatrix(t, 1) = T_outside; % Tracks outdoor temperatures

    Qleftgraph(t) = Qleft; % Tracks Q flux left side
    Qrightgraph(t) = Qright; % Tracks Q flux right side
    q_net(t) = (Qleft) - (Qright); % Tracks net Q flux through the wall


end



% Print PCM data
fprintf('\nThickness of PCM: %.4f m', Dxt);
fprintf('\nPCM Mass: %.2f kg', Mtotal);
fprintf('\nPCM energy density: %.3f kJ/kg', (NetQstore/1000)/Mtotal);
fprintf('\nThe total time was %.f minutes (%.f hours).', t_max / 60, t_max / (60*60))
fprintf('\nThe total energy stored over %.f minutes was %.4f kJ.\n ', t / 60, NetQstore / 1000);

% Display variables
dateFormatted = strrep(dateStr, '.', '/'); % Changes input to be more display friendly
avgPCMTemp = mean(PCMTemps); % Average temperature in °C 
maxRateChange = max(abs(diff(PCMTemps))) * 3600; % Max rate of change of temperature°C/hr
totalEnergyStored = NetQstore / 1000; % kJ

% Plotting
figure;
tiledlayout(3,2);
sgtitle(['Simulation for ' upper(city) ' on ' dateFormatted ' using (' upper(PCMname) ')']);

% Plot 1
nexttile;
plot((1:time_steps)/3600, MassPercent * 100, 'b', 'LineWidth', 2);
title([upper(PCMname) ' Melted (%)']);

% Plot 2
nexttile;
plot((1:time_steps)/3600, PCMTemps, 'b', 'LineWidth', 2);
title([upper(PCMname) ' Temperature']);

% Plot 3
nexttile;
plot((1:time_steps)/3600, OutsideTempMatrix, 'b', 'LineWidth', 2);
title('Outside Temperature');

% Plot 4 - Heat Flux
nexttile;
plot((1:time_steps)/3600, q_net / A, 'b', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Heat Flux (W/m²)');
title('Heat Flux Through Wall Layers');

% Summary Stats
nexttile;
axis off;
title('Summary Statistics');

statsText = {
    sprintf('Final Melted:        %.2f %%', min(100, MassPercent(end)*100))
    sprintf('Total Energy Stored: %.2f kJ', totalEnergyStored)
    sprintf('Avg PCM Temp:        %.2f °C', avgPCMTemp)
    sprintf('Max PCM Temp:        %.2f °C', max(PCMTemps))
    sprintf('Max Temp Rate:       %.2f °C/hr', maxRateChange)
};

xpos = 0;
ypos = 0.9;
dy = 0.2;

for i = 1:length(statsText)
    text(xpos, ypos - (i-1)*dy, statsText{i}, 'FontName', 'Consolas', 'FontSize', 10);
end

viewdata = table((1:time_steps)', PCMTemps', MassPercent'*100, Qleftgraph', Qrightgraph', q_net', OutsideTempMatrix, ...
    'VariableNames', {'Time_s', 'PCM_Temp_C', 'Melted_percent', 'Q_left_W', 'Q_right_W', 'Net_Q_W', 'Outside_Temp_C'});






