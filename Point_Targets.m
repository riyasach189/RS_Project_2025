% Riya Sachdeva (2022411) Sarthak Kalpasi (2021197)
% 07/05/2025
% RS Project 2025
% SAR Stripmap Simulation for Point Targets

clear; close all; clc;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radar Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Speed of light (m/s)
c_mps = 3e8;      

% Radar center frequency (Hz)
fc_Hz = 1e9;       

% Wavelength (m)
wavelength_m = c_mps / fc_Hz;     

% Platform altitude (m)
altitude_m = 1000;        

% Platform velocity (m/s)
vp_mps = 100;        

% Squint angle (degrees)
squint_angle_degrees = 0;    

% Transmit power (W)
Ptx_watts = 50e3;       

% Antenna gain (dB)
antenna_gain_dB = 40;                        

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Waveform Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pulse Repetition Frequency (Hz)
PRF_Hz = 1000;

% Pulse Repetition Interval (s)
PRI_s = 1 / PRF_Hz;

% Duty Cycle (no units)
duty_cycle_nu = 0.1;

% Pulse Width (s)
T_s = PRI_s * duty_cycle_nu;

% Range Resolution (m)
rangeres_m = 30;

% Cross Range Resolution (m)
crossrangeres_m = 30;

% Chirp rate (Hz/s)
Kfast_Hzpsec = c_mps / (2 * rangeres_m * T_s);

% Chirp bandwidth (Hz)
bandwidth_fast_Hz = Kfast_Hzpsec * T_s;

% Sampling Frequency (Hz)
Fs_Hz = 3 * bandwidth_fast_Hz;        

% Fast-time sample interval (s)
delt_s = 1 / Fs_Hz;    

% Fast-time axis (INPUT) (s)
fast_time_axis_s = 0:delt_s:PRI_s-delt_s;  

% Number of fast-time samples (no units)
Nfast_nu = length(fast_time_axis_s);          

% Fast-time axis (OUTPUT) (s)
fast_time_axis_xcorr_s = -PRI_s+delt_s:delt_s:PRI_s-delt_s;  

% Number of fast-time samples after Matched Filtering (no units)
Nrange_nu = length(fast_time_axis_xcorr_s);   

% Range axis (m)
range_axis_m = c_mps * 0.5 * fast_time_axis_xcorr_s; 

% Coherent Processing Interval (s)
CPI_s = 5 * PRI_s; 

% Slow-time axis (INPUT) (s)
slow_time_axis_s = 0:PRI_s:CPI_s-PRI_s;  

% Number of pulses (no units)
Nslow_nu = length(slow_time_axis_s);          

% Doppler frequency resolution (Hz)
delf_Hz = 1 / CPI_s;   

% Slow-time Axis (OUTPUT) (s)
slow_time_axis_xcorr_s = -CPI_s+PRI_s:PRI_s:CPI_s-PRI_s;

% Number of cross range samples (no units)
Ncross_range_nu = length(slow_time_axis_xcorr_s);

% Cross Range Axis (m)
cross_range_axis_m = vp_mps * slow_time_axis_xcorr_s; 

% Number of samples in a pulse width (no units)
Nsamp_pw_nu = round(Fs_Hz * T_s);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitted Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define LFM waveform using Radar Toolbox
tx_waveform = phased.LinearFMWaveform('SampleRate', Fs_Hz, ...
    'PulseWidth', T_s, ...
    'PRF', PRF_Hz, ...
    'SweepBandwidth', bandwidth_fast_Hz, ...
    'SweepDirection','Up');

% Transmitted signal (volts)
stx_volts = tx_waveform();                    

% Plot transmitted signal (real and imaginary parts)
figure;
plot(fast_time_axis_s * 1e3, real(stx_volts)); 
xlabel('Time (ms)'); ylabel('Amplitude');
title('Real Part of Transmitted Signal');
grid on;

figure;
plot(fast_time_axis_s * 1e3, imag(stx_volts)); 
xlabel('Time (ms)'); ylabel('Amplitude');
title('Imaginary Part of Transmitted Signal');
grid on;

% Zoomed Time Axis (small part of entire time axis) (seconds) (for Tx signal)
time_axis_zoomed_tx_s = fast_time_axis_s(1:round(Nsamp_pw_nu));

% Zoomed Transmitted Signal (zoomed at region where signal is non-zero) (volts)
stx_zoomed_volts = stx_volts(1:round(Nsamp_pw_nu));

% Transmitted Signal Plot (zoomed at region where signal is non-zero)
figure;
plot(time_axis_zoomed_tx_s, real(stx_zoomed_volts));
xlabel('Time Axis - Zoomed (seconds)');
ylabel('Transmitted Signal (volts)');
title('Real Part of Transmitted Signal (Zoomed Plot)');

figure;
plot(time_axis_zoomed_tx_s, imag(stx_zoomed_volts));
xlabel('Time Axis - Zoomed (seconds)');
ylabel('Transmitted Signal (volts)');
title('Imaginary Part of Transmitted Signal (Zoomed Plot)');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antenna and Platform Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Isotropic TX antenna
tx_antenna = phased.IsotropicAntennaElement('FrequencyRange', [1e9 20e9]);  

% Isotropic RX Antenna
rx_antenna = clone(tx_antenna);   

% Transmitter
transmitter = phased.Transmitter('PeakPower', Ptx_watts, 'Gain', antenna_gain_dB);

% Radiator
radiator = phased.Radiator('Sensor', tx_antenna,'OperatingFrequency', fc_Hz, 'PropagationSpeed', c_mps);

% Collector
collector = phased.Collector('Sensor', rx_antenna,'OperatingFrequency', fc_Hz, 'PropagationSpeed', c_mps);

% Receiver
receiver = phased.ReceiverPreamp('Gain', antenna_gain_dB, 'NoiseFigure', 30);     

% Moving Radar Platform
radar_platform = phased.Platform('InitialPosition', [0; -200; altitude_m], 'Velocity', [0; vp_mps; 0]);          

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target Parameters (3 Point targets)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free Space Propagation Channel
channel = phased.FreeSpace('SampleRate', Fs_Hz, ...
    'OperatingFrequency', fc_Hz, ...
    'TwoWayPropagation', true, ...
    'PropagationSpeed', c_mps);   

% Number of Targets (no units)
num_tgts = 3;

% Target Positions (m)
targetpos_m = [800,0,0; 1000,0,0; 1300,0,0]'; 

% Target Velocities (mps)
targetvel_mps = [0,0,0; 0,0,0; 0,0,0]';

% Point Targets with RCS = 1 m^2
target = phased.RadarTarget('MeanRCS', [1, 1, 1], 'OperatingFrequency', fc_Hz);
pointTargets = phased.Platform('InitialPosition', targetpos_m, 'Velocity', targetvel_mps);

% Plot Ground Truth of Target Positions
figure;
hold on; grid on; axis equal;

% Plot targets
plot3(targetpos_m(1, :), targetpos_m(2, :), targetpos_m(3, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
text(targetpos_m(1, :)+20, targetpos_m(2, :)+20, targetpos_m(3, :), {'Target 1','Target 2','Target 3'});

% Plot radar platform initial position
radar_init_pos = radar_platform.InitialPosition;
plot3(radar_init_pos(1), radar_init_pos(2), radar_init_pos(3), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
text(radar_init_pos(1)+20, radar_init_pos(2)+20, radar_init_pos(3), 'Radar');

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Ground Truth: Radar and Target Positions');
legend('Targets', 'Radar Platform');
view(3);  % 3D view
hold off;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Received Signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
srx_volts = zeros(Nfast_nu, Nslow_nu);

for n = 1:Nslow_nu
    % Update radar platform and target positions
    [radar_pos, radar_vel] = radar_platform(PRI_s);
    [tgt_pos, tgt_vel] = pointTargets(PRI_s);

    % Get the range and angle to the point targets
    [~, tgt_angle] = rangeangle(tgt_pos, radar_pos);

    % Generate the waveform and truncate to useful length
    txsig = tx_waveform();
    txsig = txsig(1:Nfast_nu); 

    % Transmit signal
    txsig = transmitter(txsig);

    % Set azimuth to 0 for broadside
    tgt_angle(1,:) = 0;

    % Radiate signal
    txsig = radiator(txsig, tgt_angle);

    % Propagate to targets
    txsig = channel(txsig, radar_pos, tgt_pos, radar_vel, tgt_vel);

    % Reflect from targets 
    txsig = target(txsig);

    % Collect return
    txsig = collector(txsig, tgt_angle);

    % Receive signal
    srx_volts(:, n) = receiver(txsig);
end

% Plot received signal 
figure;
imagesc(slow_time_axis_s, fast_time_axis_s, real(srx_volts));
xlabel('Slow Time (s)');
ylabel('Fast Time (s)');
title('Real Part of Received Signal');
colorbar;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matched Filtering (Pulse Compression)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MF_Mat = zeros(Nrange_nu, Nslow_nu);

for n = 1:Nslow_nu
    MF_Mat(:, n) = xcorr(srx_volts(:, n), stx_volts);
end

% Plot SAR range-compressed data
figure;
imagesc(slow_time_axis_s, range_axis_m, abs(MF_Mat));
xlabel('Slow Time (s)');
ylabel('Range (m)');
title('SAR Range Compressed Data');
colormap('jet');
colorbar;
ylim([0 2500]);

