% Riya Sachdeva (2022411) Sarthak Kalpasi (2021197)
% 05/05/2025
% RS Project 2025
% SAR Stripmap Simulation for Point Targets and Terrain

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

% Initialize a matrix to store the collected raw IQ data
raw = zeros(rangeSamples, numPulses);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Collection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Collecting radar data (%d pulses)...\n', numPulses);

% Initialize a counter for collected pulses
ii = 1;

% Create a figure to visualize the raw IQ data as it is collected
hRaw = helperPlotRawIQ(raw,minSample);

% Advance the scenario and collect data for each pulse
while advance(scene)
    % Receive radar data for the current time step
    tmp = receive(scene);
    % Extract the relevant range samples and store in the raw data matrix
    raw(:,ii) = tmp{1}(minSample:truncRngSamp);

    % Update the plot periodically to show data collection progress
    if mod(ii,5) == 0
        helperUpdatePlotRawIQ(hRaw,raw);
    end

    % Increment the pulse counter
    ii = ii + 1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAR Image Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Processing SAR image...\n');

% Perform range migration processing on the raw IQ data to create an SLC image
slcimg = rangeMigrationLFM(raw,rdr.Waveform,freq_Hz,v_mps,rc_m);

% Plot the generated SLC image alongside the ground truth
helperPlotSLC(slcimg,minSample,fs_Hz,v_mps,prf_Hz,rdrpos1,targetpos,xvec,yvec,A);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target Visibility Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nTarget Visibility Summary:\n');

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
    figure('Position',[750 100 600 400]);
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
