% Riya Sachdeva (2022411) Sarthak Kalpasi (2021197)
% 05/07/2025
% RS Project 2025
% SAR Stripmap Simulation for Point Targets and Terrain with Weather Effects

clear; close all; clc;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Terrain Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize random number generator for reproducibility
rng(2021);

% Define limits for the terrain in X and Y coordinates
xLimits = [900 1200];
yLimits = [-200 200];

% Set parameters for terrain generation using the Diamond-Square algorithm
roughnessFactor = 3;      % Controls the "bumpiness" of the terrain
initialHeight = 0;        % Starting height for the terrain
initialPerturb = 200;     % Initial random perturbation
numIter = 4;              % Number of iterations for refinement

% Generate the terrain data
fprintf('Generating terrain...\n');
[x,y,A] = helperRandomTerrainGenerator(roughnessFactor,initialHeight,initialPerturb,...
    xLimits(1),xLimits(2),yLimits(1),yLimits(2),numIter);

% Ensure terrain height is non-negative
A(A < 0) = 0;

% Extract X and Y vectors from the meshgrid
xvec = x(1,:);
yvec = y(:,1);

% Plot the simulated terrain for visualization
helperPlotSimulatedTerrain(xvec,yvec,A);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radar Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define radar center frequency (Hz)
freq_Hz = 1e9;

% Calculate wavelength (m) and speed of light (m/s) from frequency
[lambda_m,c_mps] = freq2wavelen(freq_Hz);

% Define radar bandwidth (Hz)
bw_Hz = 30e6;

% Define radar sampling frequency (Hz)
fs_Hz = 60e6;

% Define transmit pulse duration (s)
tpd_s = 3e-6;

% Calculate range resolution (m) from bandwidth
rngRes_m = bw2rangeres(bw_Hz);

% Define radar aperture length (m)
apertureLength_m = 6;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Platform Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimized platform velocity (m/s)
v_mps = 100;

% Data collection duration (s)
dur_s = 0.3;

% Radar platform altitude (m)
rdrhgt_m = 1000;

% Initial radar platform position (m)
rdrpos1 = [0 0 rdrhgt_m];

% Radar platform velocity vector (m/s)
rdrvel = [0 v_mps 0];

% Final radar platform position after the data collection duration (m)
rdrpos2 = rdrvel*dur_s + rdrpos1;

% Calculate SAR aperture length (m)
len_m = sarlen(v_mps,dur_s);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define target positions (m)
targetpos = [1000,len_m/2,0;1020,len_m/2,0;1040,len_m/2,0];

% Initialize target heights (m)
tgthgts = 110*ones(1,3);

% Adjust target heights based on terrain elevation at their locations
for it = 1:3
    % Find the nearest terrain point for the target's X coordinate
    [~,idxX] = min(abs(targetpos(it,1) - xvec));
    % Find the nearest terrain point for the target's Y coordinate
    [~,idxY] = min(abs(targetpos(it,2) - yvec));
    % Add terrain height to the base target height
    tgthgts(it) = tgthgts(it) + A(idxX,idxY);
    % Update the target's Z coordinate (height)
    targetpos(it,3) = tgthgts(it);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAR Processing Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the center range (m) to the targets
rc_m = sqrt((rdrhgt_m - mean(tgthgts))^2 + (mean(targetpos(:,1)))^2);

% Calculate the depression angle (degrees) to the targets
depang_deg = depressionang(rdrhgt_m,rc_m,'Flat','TargetHeight',mean(tgthgts));

% Grazing angle is approximately equal to the depression angle in this scenario
grazang_deg = depang_deg;

% Calculate potential swath width and minimum/maximum range (not used later in this script but informative)
[~,~] = aperture2swath(rc_m,lambda_m,apertureLength_m,grazang_deg);

% Calculate PRF bounds (not used later in this script but informative)
[~,~] = sarprfbounds(v_mps,bw2rangeres(bw_Hz),len_m,grazang_deg);

% Set the Pulse Repetition Frequency (PRF) to create an integer ratio with fs
% This simplifies some processing steps
prf_Hz = 300;

% Display the chosen PRF and check if the ratio with fs is an integer
fprintf('PRF: %d Hz, fs/PRF ratio: %d (integer check: %d)\n', ...
    prf_Hz, fs_Hz/prf_Hz, mod(fs_Hz/prf_Hz,1)==0);

% Plot the ground truth: terrain, radar path, and target positions
helperPlotGroundTruth(xvec,yvec,A,rdrpos1,rdrpos2,targetpos);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reflectivity Map Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define grazing angle and frequency tables for reflectivity data
grazTable = 20:1:60;
freqTable = [1e9 10e9];

% Define the number of different surface types
numSurfaces = 2;

% Initialize a matrix to store reflectivity data for each surface type
reflectivityLayers = zeros(numel(grazTable),numel(freqTable),numSurfaces);

% Load land reflectivity data for "Woods"
reflectivityLayers(:,:,1) = landreflectivity('Woods', grazTable,freqTable);

% Load land reflectivity data for "WoodedHills"
reflectivityLayers(:,:,2) = landreflectivity('WoodedHills', grazTable,freqTable);

% Create a reflectivity type map based on terrain height
% Terrain above 100m is considered "WoodedHills" (type 2), otherwise "Woods" (type 1)
reflectivityType = ones(size(A));
reflectivityType(A > 100) = 2;

% Plot the reflectivity map over the terrain
helperPlotReflectivityMap(xvec,yvec,A,reflectivityType,rdrpos1,rdrpos2,targetpos);

% Create a surfaceReflectivity object for the radar scenario
reflectivityMap = surfaceReflectivity('Custom','Frequency',freqTable,...
    'GrazingAngle',grazTable,'Reflectivity',reflectivityLayers,...
    'Speckle','Rayleigh'); % Include Rayleigh speckle model

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weather Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define weather condition intensities
% Rain rates in mm/hour
rainRates = [1, 200, 1000];  % Low, Medium, High
% Snow rates in mm/hour
snowRates = [5, 50, 300];  % Low, Medium, High
% Fog visibility in meters
fogVisibility = [500, 200, 50];  % Low, Medium, High

% Names for labeling
weatherNames = {'Clear', 'Rain-Low', 'Rain-Medium', 'Rain-High', ...
                'Snow-Low', 'Snow-Medium', 'Snow-High', ...
                'Fog-Low', 'Fog-Medium', 'Fog-High'};

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radar Scenario Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Creating radar scenario...\n');

% Create a radar scenario object
scene = radarScenario('UpdateRate',prf_Hz,'IsEarthCentered',false,'StopTime',dur_s);

% Add the radar platform to the scenario
rdrplat = platform(scene,'Trajectory',kinematicTrajectory('Position',rdrpos1,'Velocity',[0 v_mps 0]));

% Add the point targets to the scenario
rcs = rcsSignature('Pattern',5); % Define an RCS signature (constant in this case)
for it = 1:3
    platform(scene,'Position',targetpos(it,:),'Signatures',{rcs});
end

% Add the land surface to the scenario, including terrain and reflectivity
s = landSurface(scene,'Terrain',A,'Boundary',[xLimits;yLimits],...
    'RadarReflectivity',reflectivityMap,...
    'ReflectivityMap',reflectivityType);

% Create a radar transceiver object
maxRange_m = 2500; % Maximum range for radar operation
mountAngles = [0 depang_deg 0]; % Mounting angles [azimuth, elevation, roll]
rdr = radarTransceiver('MountingAngles',mountAngles,'NumRepetitions',1,...
    'RangeLimits',[0 maxRange_m]);

% Set radar transceiver parameters
rdr.Transmitter.PeakPower = 50e3; % Peak transmit power (W)
rdr.Receiver.SampleRate = fs_Hz; % Receiver sample rate (Hz)
rdr.Receiver.NoiseFigure = 30; % Receiver noise figure (dB)

% Define the radar antenna using a Sinc antenna element
antbw_deg = ap2beamwidth(apertureLength_m,lambda_m); % Calculate antenna beamwidth
ant = phased.SincAntennaElement('FrequencyRange',[1e9 10e9],'Beamwidth',antbw_deg);

% Assign the antenna to the transmit and receive paths
rdr.TransmitAntenna.Sensor = ant;
rdr.TransmitAntenna.OperatingFrequency = freq_Hz;
rdr.ReceiveAntenna.Sensor = ant;
rdr.ReceiveAntenna.OperatingFrequency = freq_Hz;

% Calculate and set antenna gain
antennaGain = aperture2gain(apertureLength_m^2,lambda_m);
rdr.Transmitter.Gain = antennaGain;
rdr.Receiver.Gain = antennaGain;

% Configure the Linear FM waveform for the radar
rdr.Waveform = phased.LinearFMWaveform('SampleRate',fs_Hz,'PulseWidth',tpd_s,...
    'PRF',prf_Hz,'SweepBandwidth',bw_Hz);

% Add the radar transceiver to the radar platform
rdrplat.Sensors = rdr;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clutter Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Generating clutter...\n');

% Configure clutter generation within the scenario
clutterGenerator(scene,rdr,'Resolution',rngRes_m,'RangeLimit',maxRange_m);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IQ Data Collection Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the minimum sample for collecting data (to exclude direct blast)
minSample = 500;

% Calculate the maximum range time and corresponding sample index
maxRange_time = range2time(maxRange_m);
truncRngSamp = ceil(maxRange_time*fs_Hz);

% Calculate the number of range samples to collect
rangeSamples = floor(truncRngSamp - minSample + 1);

% Calculate the pulse repetition interval (s)
T = 1/prf_Hz;

% Calculate the total number of pulses to collect
numPulses = floor(dur_s/T + 1);

% Display the number of range samples and pulses
fprintf('Range samples: %d, Number of pulses: %d\n', rangeSamples, numPulses);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Collection & Processing for All Weather Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize a cell array to store raw IQ data for each weather condition
rawData = cell(1, 10); % 1 for clear + 9 for weather conditions
slcImages = cell(1, 10); % To store processed SLC images

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Weather Data Collection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Collecting radar data (clear weather)...\n');
rawData{1} = zeros(rangeSamples, numPulses);
ii = 1;

% Create a figure for raw IQ data visualization
figure('Position',[750 100 600 400]);
hRawPlot = helperPlotRawIQ(rawData{1}, minSample);

while advance(scene)
    tmp = receive(scene);
    rawData{1}(:,ii) = tmp{1}(minSample:truncRngSamp);
    
    if mod(ii,5) == 0
        helperUpdatePlotRawIQ(hRawPlot, rawData{1});
    end
    
    ii = ii + 1;
end

% Process clear weather data
slcImages{1} = rangeMigrationLFM(rawData{1}, rdr.Waveform, freq_Hz, v_mps, rc_m);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weather Effect Data Collection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate approximate range to targets (for weather attenuation)
targetRange = mean(sqrt(sum((targetpos - rdrpos1).^2, 2)));

% Rain (Low, Medium, High)
for i = 1:3
    weatherIdx = i + 1;
    fprintf('Collecting data for Rain (%d mm/hr)...\n', rainRates(i));
    
    % Create a new scenario for each weather condition
    weatherScene = radarScenario('UpdateRate', prf_Hz, 'IsEarthCentered', false, 'StopTime', dur_s);
    
    % Add platform, targets and land surface as before
    weatherPlatform = platform(weatherScene, 'Trajectory', kinematicTrajectory('Position', rdrpos1, 'Velocity', [0 v_mps 0]));
    weatherPlatform.Sensors = rdr;
    
    for it = 1:3
        platform(weatherScene, 'Position', targetpos(it,:), 'Signatures', {rcs});
    end
    
    landSurface(weatherScene, 'Terrain', A, 'Boundary', [xLimits; yLimits], ...
        'RadarReflectivity', reflectivityMap, 'ReflectivityMap', reflectivityType);
    
    % Configure clutter for the new scenario
    clutterGenerator(weatherScene, rdr, 'Resolution', rngRes_m, 'RangeLimit', maxRange_m);
    
    % Collect data for this weather condition
    rawData{weatherIdx} = zeros(rangeSamples, numPulses);
    pulseIdx = 1;
    
    while advance(weatherScene)
        % Get received signal
        tmp = receive(weatherScene);
        sampleData = tmp{1}(minSample:truncRngSamp);
        
        % Apply rain attenuation - using a simplified model instead of rainpl
        % Simplified attenuation model: 0.01 dB/km per mm/hr at 1 GHz
        attenuationDb = 0.01 * rainRates(i) * (targetRange/1000);
        attenuationFactor = 10^(-attenuationDb/20);
        rawData{weatherIdx}(:, pulseIdx) = sampleData * attenuationFactor;
        
        pulseIdx = pulseIdx + 1;
    end
    
    % Process the weather-affected data
    slcImages{weatherIdx} = rangeMigrationLFM(rawData{weatherIdx}, rdr.Waveform, freq_Hz, v_mps, rc_m);
end

% Snow (Low, Medium, High)
for i = 1:3
    weatherIdx = i + 4;
    fprintf('Collecting data for Snow (%g mm/hr)...\n', snowRates(i));
    
    % Create a new scenario for each weather condition
    weatherScene = radarScenario('UpdateRate', prf_Hz, 'IsEarthCentered', false, 'StopTime', dur_s);
    
    % Add platform, targets and land surface as before
    weatherPlatform = platform(weatherScene, 'Trajectory', kinematicTrajectory('Position', rdrpos1, 'Velocity', [0 v_mps 0]));
    weatherPlatform.Sensors = rdr;
    
    for it = 1:3
        platform(weatherScene, 'Position', targetpos(it,:), 'Signatures', {rcs});
    end
    
    landSurface(weatherScene, 'Terrain', A, 'Boundary', [xLimits; yLimits], ...
        'RadarReflectivity', reflectivityMap, 'ReflectivityMap', reflectivityType);
    
    % Configure clutter for the new scenario
    clutterGenerator(weatherScene, rdr, 'Resolution', rngRes_m, 'RangeLimit', maxRange_m);
    
    % Collect data for this weather condition
    rawData{weatherIdx} = zeros(rangeSamples, numPulses);
    pulseIdx = 1;
    
    while advance(weatherScene)
        % Get received signal
        tmp = receive(weatherScene);
        sampleData = tmp{1}(minSample:truncRngSamp);
        
        % Apply snow attenuation - using a simplified model
        % Simplified attenuation: 0.05 dB/km per mm/hr for snow at 1 GHz
        attenuationDb = 0.05 * snowRates(i) * (targetRange/1000);
        attenuationFactor = 10^(-attenuationDb/20);
        rawData{weatherIdx}(:, pulseIdx) = sampleData * attenuationFactor;
        
        pulseIdx = pulseIdx + 1;
    end
    
    % Process the weather-affected data
    slcImages{weatherIdx} = rangeMigrationLFM(rawData{weatherIdx}, rdr.Waveform, freq_Hz, v_mps, rc_m);
end

% Fog (Low, Medium, High)
for i = 1:3
    weatherIdx = i + 7;
    fprintf('Collecting data for Fog (%d m visibility)...\n', fogVisibility(i));
    
    % Create a new scenario for each weather condition
    weatherScene = radarScenario('UpdateRate', prf_Hz, 'IsEarthCentered', false, 'StopTime', dur_s);
    
    % Add platform, targets and land surface as before
    weatherPlatform = platform(weatherScene, 'Trajectory', kinematicTrajectory('Position', rdrpos1, 'Velocity', [0 v_mps 0]));
    weatherPlatform.Sensors = rdr;
    
    for it = 1:3
        platform(weatherScene, 'Position', targetpos(it,:), 'Signatures', {rcs});
    end
    
    landSurface(weatherScene, 'Terrain', A, 'Boundary', [xLimits; yLimits], ...
        'RadarReflectivity', reflectivityMap, 'ReflectivityMap', reflectivityType);
    
    % Configure clutter for the new scenario
    clutterGenerator(weatherScene, rdr, 'Resolution', rngRes_m, 'RangeLimit', maxRange_m);
    
    % Collect data for this weather condition
    rawData{weatherIdx} = zeros(rangeSamples, numPulses);
    pulseIdx = 1;
    
    while advance(weatherScene)
        % Get received signal
        tmp = receive(weatherScene);
        sampleData = tmp{1}(minSample:truncRngSamp);
        
        % Apply fog attenuation - using a simplified model
        % Simplified model: 0.4 * freq_GHz^2 / visibility_km dB/km
        visibility_km = fogVisibility(i)/1000;  % Convert m to km
        attenuationDb = 0.4 * (freq_Hz/1e9)^2 / visibility_km * (targetRange/1000);
        attenuationFactor = 10^(-attenuationDb/20);
        rawData{weatherIdx}(:, pulseIdx) = sampleData * attenuationFactor;
        
        pulseIdx = pulseIdx + 1;
    end
    
    % Process the weather-affected data
    slcImages{weatherIdx} = rangeMigrationLFM(rawData{weatherIdx}, rdr.Waveform, freq_Hz, v_mps, rc_m);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot SAR Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot SAR image for clear weather (original)
helperPlotSLC(slcImages{1}, minSample, fs_Hz, v_mps, prf_Hz, rdrpos1, targetpos, xvec, yvec, A);

% Plot weather comparison figures
% Rain comparison
figure('Position', [100 100 1200 400]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot rain images
c = physconst('LightSpeed');
numSamples = size(slcImages{1}, 1);
samples = minSample:(numSamples + minSample - 1);
sampleTime = samples*1/fs_Hz;
rngVec = time2range(sampleTime(1:end), c);

numPulses = size(slcImages{1}, 2);
y = linspace(0, 100, numPulses);

% Find common max value for consistent colormapping
maxVal = 0;
for i = 2:4
    maxVal = max(maxVal, max(abs(slcImages{i}(:))));
end

% Rain plots
for i = 2:4
    nexttile;
    slcimg = abs(slcImages{i}).';
    pcolor(rngVec, y, slcimg);
    shading flat;
    colormap(parula);
    clim([0 maxVal]);
    
    if i == 4
        colorbar('eastoutside');
    end
    
    xlabel('Slant Range (m)');
    if i == 2
        ylabel('Cross-range (m)');
    end
    
    title(sprintf('Rain: %d mm/hr', rainRates(i-1)));
    axis equal;
    xlim([1250 1420]);
    ylim([0 100]);
end
sgtitle('SAR Image Comparison - Rain Conditions', 'FontSize', 14);

% Snow comparison
figure('Position', [100 100 1200 400]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Find common max value for consistent colormapping
maxVal = 0;
for i = 5:7
    maxVal = max(maxVal, max(abs(slcImages{i}(:))));
end

% Snow plots
for i = 5:7
    nexttile;
    slcimg = abs(slcImages{i}).';
    pcolor(rngVec, y, slcimg);
    shading flat;
    colormap(parula);
    clim([0 maxVal]);
    
    if i == 7
        colorbar('eastoutside');
    end
    
    xlabel('Slant Range (m)');
    if i == 5
        ylabel('Cross-range (m)');
    end
    
    title(sprintf('Snow: %g mm/hr', snowRates(i-4)));
    axis equal;
    xlim([1250 1420]);
    ylim([0 100]);
end
sgtitle('SAR Image Comparison - Snow Conditions', 'FontSize', 14);

% Fog comparison
figure('Position', [100 100 1200 400]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Find common max value for consistent colormapping
maxVal = 0;
for i = 8:10
    maxVal = max(maxVal, max(abs(slcImages{i}(:))));
end

% Fog plots
for i = 8:10
    nexttile;
    slcimg = abs(slcImages{i}).';
    pcolor(rngVec, y, slcimg);
    shading flat;
    colormap(parula);
    clim([0 maxVal]);
    
    if i == 10
        colorbar('eastoutside');
    end
    
    xlabel('Slant Range (m)');
    if i == 8
        ylabel('Cross-range (m)');
    end
    
    title(sprintf('Fog: %d m visibility', fogVisibility(i-7)));
    axis equal;
    xlim([1250 1420]);
    ylim([0 100]);
end
sgtitle('SAR Image Comparison - Fog Conditions', 'FontSize', 14);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target Visibility Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nTarget Visibility Summary (Clear Weather):\n');

% Analyze the visibility of each target throughout the data collection
for it = 1:3
    % Initialize an array to track occlusion for each pulse
    occ = false(1,numPulses);
    % Iterate through each pulse to check for occlusion
    for ip = 1:numPulses
        % Calculate the radar position at the time of this pulse
        rdrpos = rdrpos1 + rdrvel.*1/prf_Hz*(ip - 1);
        % Check if the target is occluded by the terrain at this radar position
        occ(ip) = s.occlusion(rdrpos,targetpos(it,:));
    end
    % Report the target visibility status
    helperGetVisibilityStatus(it,occ);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target Visibility Analysis for All Weather Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nTarget Visibility Analysis for All Weather Conditions:\n');

% Map indices to condition names for reporting
weatherNames = {'Clear', 'Rain-Low', 'Rain-Medium', 'Rain-High', ...
                'Snow-Low', 'Snow-Medium', 'Snow-High', ...
                'Fog-Low', 'Fog-Medium', 'Fog-High'};

% Calculate range vector for target localization
c = physconst('LightSpeed');
numSamples = size(slcImages{1}, 1);
samples = minSample:(numSamples + minSample - 1);
sampleTime = samples*1/fs_Hz;
rngVec = time2range(sampleTime(1:end), c);

% Calculate cross-range vector
numPulses = size(slcImages{1}, 2);
crossRngVec = linspace(0, 100, numPulses);

% Target slant ranges (approximate)
targetSlantRanges = sqrt(targetpos(:,1).^2 + targetpos(:,3).^2);

% Find target indices in SLC images
targetIndices = zeros(3, 2); % 3 targets, 2 dimensions (range, cross-range)
for i = 1:3
    % Find closest range bin
    [~, targetIndices(i,1)] = min(abs(rngVec - targetSlantRanges(i)));
    
    % Find closest cross-range bin (assuming targets are centered in crossrange)
    [~, targetIndices(i,2)] = min(abs(crossRngVec - targetpos(i,2)));
end

% Window size for target analysis (pixels)
windowSize = 5;

% Initialize visibility results table
visibilityResults = zeros(3, length(weatherNames));

% First, get peak values for clear weather (reference)
clearImg = abs(slcImages{1});
clearPeaks = zeros(3, 1);

for t = 1:3
    % Target center indices
    tgtR = targetIndices(t,1);
    tgtC = targetIndices(t,2);
    
    % Define window around target
    rWin = max(1, tgtR-windowSize):min(size(clearImg,1), tgtR+windowSize);
    cWin = max(1, tgtC-windowSize):min(size(clearImg,2), tgtC+windowSize);
    
    % Extract target window and get peak
    targetWindow = clearImg(rWin, cWin);
    clearPeaks(t) = max(targetWindow(:));
    
    % Clear weather is 100% visible by definition
    visibilityResults(t, 1) = 100;
end

% Minimum detectability threshold (% of clear weather peak)
minDetectable = 5; 

% Analyze each weather condition
for w = 2:length(weatherNames)  % Skip clear weather (already set to 100%)
    % Get the corresponding SLC image
    slcImg = abs(slcImages{w});
    
    % Analyze each target
    for t = 1:3
        % Target center indices
        tgtR = targetIndices(t,1);
        tgtC = targetIndices(t,2);
        
        % Define window around target
        rWin = max(1, tgtR-windowSize):min(size(slcImg,1), tgtR+windowSize);
        cWin = max(1, tgtC-windowSize):min(size(slcImg,2), tgtC+windowSize);
        
        % Extract target window
        targetWindow = slcImg(rWin, cWin);
        
        % Find peak in window
        peakValue = max(targetWindow(:));
        
        % Calculate relative visibility as percentage of clear weather peak
        relativeVisibility = (peakValue / clearPeaks(t)) * 100;
        
        % If below minimum detectability, set to 0
        if relativeVisibility < minDetectable
            relativeVisibility = 0;
        end
        
        % Store result
        visibilityResults(t, w) = relativeVisibility;
    end
end

% Print visibility results
fprintf('\nTarget Visibility Results (Percentage relative to clear weather):\n');
fprintf('%-12s', 'Weather:');
for w = 1:length(weatherNames)
    fprintf('%-12s', weatherNames{w});
end
fprintf('\n');

for t = 1:3
    fprintf('Target %d:   ', t);
    for w = 1:length(weatherNames)
        fprintf('%-10.1f%%', visibilityResults(t, w));
    end
    fprintf('\n');
end

% Convert to average visibility by condition
fprintf('\nAverage Target Visibility By Weather Condition:\n');
for w = 1:length(weatherNames)
    avgVisibility = mean(visibilityResults(:, w));
    fprintf('%s: %.1f%% visibility compared to clear weather\n', weatherNames{w}, avgVisibility);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions (Provided)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Helper function to generate random terrain ---
function [x,y,terrain] = helperRandomTerrainGenerator(f,initialHeight,initialPerturb,minX,maxX,minY,maxY,numIter)
    dX = (maxX-minX)/2;
    dY = (maxY-minY)/2;
    [x,y] = meshgrid(minX:dX:maxX,minY:dY:maxY);
    terrain = ones(3,3)*initialHeight;
    perturb = initialPerturb;
    for ii = 2:numIter
        perturb = perturb/f;
        oldX = x;
        oldY = y;
        dX = (maxX-minX)/2^ii;
        dY = (maxY-minY)/2^ii;
        [x,y] = meshgrid(minX:dX:maxX,minY:dY:maxY);
        terrain = griddata(oldX,oldY,terrain,x,y);
        terrain = terrain + perturb*randn(1+2^ii,1+2^ii);
        terrain(terrain < 0) = 0;
    end
end

% --- Helper function to create a colormap for land elevation ---
function cmap = landColormap(n)
    c = hsv2rgb([5/12 1 0.4; 0.25 0.2 1; 5/72 1 0.4]);
    cmap = zeros(n,3);
    cmap(:,1) = interp1(1:3,c(:,1),linspace(1,3,n));
    cmap(:,2) = interp1(1:3,c(:,2),linspace(1,3,n));
    cmap(:,3) = interp1(1:3,c(:,3),linspace(1,3,n));
end

% --- Helper function to plot the simulated terrain ---
function helperPlotSimulatedTerrain(xvec,yvec,A)
    figure('Position',[100 550 600 400]);
    hS = surf(xvec,yvec,A);
    hS.EdgeColor = 'none';
    hC = colorbar;
    hC.Label.String = 'Elevation (m)';
    colormap(landColormap(64));
    xlabel('X (m)')
    ylabel('Y (m)')
    axis equal;
    title('Simulated Terrain')
    view([78 78])
end

% --- Helper function to plot the ground truth scene ---
function helperPlotGroundTruth(xvec,yvec,A,rdrpos1,rdrpos2,targetpos)
    figure('Position',[100 100 600 400]);
    hLim = surf([0 1200],[-200 200].',-100*ones(2),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.7);
    hold on;
    hS = surf(xvec,yvec,A);
    hS.EdgeColor = 'none';
    hC = colorbar;
    hC.Label.String = 'Elevation (m)';
    colormap(landColormap(64));
    hPlatPath = plot3([rdrpos1(1) rdrpos2(1)],[rdrpos1(2) rdrpos2(2)],[rdrpos1(3) rdrpos2(3)], ...
        '-k','LineWidth',2);
    hPlatStart = plot3(rdrpos1(1),rdrpos1(2),rdrpos1(3), ...
        'o','LineWidth',2,'MarkerFaceColor','g','MarkerEdgeColor','k');
    hTgt = plot3(targetpos(:,1),targetpos(:,2),targetpos(:,3), ...
        'o','LineWidth',2,'MarkerFaceColor',[0.8500 0.3250 0.0980], ...
        'MarkerEdgeColor','k');
    view([26 75])
    xlabel('Range (m)')
    ylabel('Cross-range (m)')
    title('Ground Truth')
    axis tight;
    zlim([-100 1200])
    legend([hLim,hS,hPlatPath,hPlatStart,hTgt], ...
        {'Scene Limits','Terrain','Radar Path','Radar Start','Target'},'Location','SouthWest')
end

% --- Helper function to plot the reflectivity map ---
function helperPlotReflectivityMap(xvec,yvec,A,reflectivityType,rdrpos1,rdrpos2,targetpos)
    figure('Position',[750 550 600 400]);
    hLim = surf([0 1200],[-200 200].',-100*ones(2),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.7);
    hold on
    hS = surf(xvec,yvec,A,reflectivityType);
    hS.EdgeColor = 'none';
    hold on;
    colormap(summer(2));
    hC = colorbar;
    clim([1 2]);
    hC.Ticks = [1 2];
    hC.TickLabels = {'Woods','Hills'};
    hC.Label.String = 'Land Type';
    hPlatPath = plot3([rdrpos1(1) rdrpos2(1)],[rdrpos1(2) rdrpos2(2)],[rdrpos1(3) rdrpos2(3)], ...
        '-k','LineWidth',2);
    hPlatStart = plot3(rdrpos1(1),rdrpos1(2),rdrpos1(3), ...
        'o','MarkerFaceColor','g','MarkerEdgeColor','k');
    hTgt = plot3(targetpos(:,1),targetpos(:,2),targetpos(:,3), ...
        'o','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor','k');
    view([26 75])
    xlabel('X (m)')
    ylabel('Y (m)')
    title('Reflectivity Map')
    axis tight;
    zlim([-100 1200])
    legend([hLim,hS,hPlatPath,hPlatStart,hTgt], ...
        {'Scene Limits','Reflectivity Map','Radar Path','Radar Start','Target'},'Location','SouthWest')
end

% --- Helper function to create a plot for raw IQ data visualization ---
function hRaw = helperPlotRawIQ(raw,minSample)
    [m,n] = size(raw);
    hRaw = pcolor(minSample:(m + minSample - 1),1:n,real(raw.'));
    hRaw.EdgeColor = 'none';
    title('Raw Data')
    xlabel('Range Samples')
    ylabel('Cross-range Samples')
    hC = colorbar;
    clim([-0.06 0.06])
    hC.Label.String = 'real(IQ)';
end

% --- Helper function to update the raw IQ data visualization plot ---
function helperUpdatePlotRawIQ(hRaw,raw)
    hRaw.CData = real(raw.');
    clim([-0.06 0.06]);
    drawnow limitrate;
end

% --- Helper function to plot the SLC image and ground truth side-by-side ---
function helperPlotSLC(slcimg,minSample,fs,v,~,~,targetpos,xvec,yvec,A)
    figure('Position',[250 300 1000 500]);
    tiledlayout(1,2,'TileSpacing','Compact');

    % Plot Ground Truth
    nexttile;
    hS = surf(xvec,yvec,A);
    hS.EdgeColor = 'none';
    hold on;
    plot3(targetpos(:,1),targetpos(:,2),targetpos(:,3), ...
        'o','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor','k');
    colormap(landColormap(64));
    hC = colorbar('southoutside');
    hC.Label.String = 'Elevation (m)';
    view([-1 75])
    xlabel('Range (m)')
    ylabel('Cross-range (m)')
    title('Ground Truth')
    axis equal
    xlim([950 1100])
    ylim([0 100])

    % Plot SAR Image
    nexttile;
    % Range vector
    c = physconst('LightSpeed');
    numSamples = size(slcimg,1);
    samples = minSample:(numSamples + minSample - 1);
    sampleTime = samples*1/fs;
    rngVec = time2range(sampleTime(1:end),c);

    % Cross-range vector - map to full 0-100m range
    numPulses = size(slcimg,2);
    y = linspace(0, 100, numPulses);

    slcimg = abs(slcimg).';
    hProc = pcolor(rngVec,y,slcimg);
    hProc.EdgeColor = 'none';
    colormap(hProc.Parent,parula)
    hC = colorbar('southoutside');
    hC.Label.String = 'Magnitude';
    xlabel('Slant Range (m)')
    ylabel('Cross-range (m)')
    title('SAR Image')
    axis equal
    xlim([1250 1420])
    ylim([0 100])
end

% --- Helper function to get and display target visibility status ---
function helperGetVisibilityStatus(tgtNum,occ)
    visibility = {'not','partially','fully'};
    if all(occ)
        idx = 1; % All pulses occluded
    elseif any(occ)
        idx = 2; % Some pulses occluded
    else
        idx = 3; % No pulses occluded
    end
    visString = visibility{idx};
    pctCollect = sum(double(~occ))./numel(occ)*100; % Percentage of time visible
    fprintf('Target %d is %s visible during the scenario (visible %.0f%% of the data collection).\n', ...
        tgtNum,visString,pctCollect)
end
