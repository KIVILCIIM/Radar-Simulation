clear all
%close all
clc

%% Radar Desgin Specifications (Skyguard Surveillance Radar)
Pd = 0.9; % Probability of detection, 99.996 denenebilir
Pfa = 1e-6; % Probability of false alarm
Max_Range= 16000; % Maximum unambiguous range of survaillance radar
RangeResolution = 150; % Required range resolution is 150m
TargetRCS = 1.0; % Required target radar cross section

%% Waveform Properties

PropagationSpeed = physconst('LightSpeed'); % Propagation speed is defined as the light speed
Pulse_BandWidth = PropagationSpeed/(2*RangeResolution); % Pulse bandwidth 
PulseWidth = 1/Pulse_BandWidth; % Pulse width
%PRF = PropagationSpeed/(2*Max_Range); % Pulse repetition frequency, it is given in the paper
PRF = 6400;
%fs = 2*Pulse_BandWidth; % Sampling rate, Nyquist–Shannon sampling theorem
fs = 6400 * 313; % fs/PRF oranının tam sayı olması gerekmektedir
waveform = phased.RectangularWaveform(...
    'PulseWidth',1/Pulse_BandWidth,...
    'PRF',PRF,...
    'SampleRate',fs);

%% Receiver Noise Characteristics
Noise_BandWidth = Pulse_BandWidth;

receiver = phased.ReceiverPreamp(...
    'Gain',40,...  %Gain yüksek seçildi
    'NoiseFigure',0,...
    'SampleRate',fs,...
    'EnableInputPort',true);

receiver.SeedSource = 'Property'; % For repeatable simulation
receiver.Seed = 2007;

%% Transmitter Characterics
num_pulse_int = 43; %survaillance radar için belirtilen sayı
SNR_Min = albersheim(Pd, Pfa, num_pulse_int);
fprintf('43 pulse için gerekli minimum SNR: %.2f dB\n', SNR_Min);
%SNR_Min=3;
Tx_Gain = 30;

fc = 10e9; % 9 GHz — typical for X-band surveillance radar
lambda = PropagationSpeed/fc;
Peak_Power = radareqpow(lambda, Max_Range, SNR_Min, PulseWidth, 'RCS', TargetRCS, 'Gain', Tx_Gain);
transmitter = phased.Transmitter('Gain',Tx_Gain, 'PeakPower',Peak_Power, 'InUseOutputPort',true);

figure;
rocsnr([-5 0 3 5 10 15], 'SignalType', 'NonfluctuatingCoherent', 'NumPulses', num_pulse_int); %ekstra değerler eklendi rocsnr içine

num_pulse_int = 43;
Tx_Gain = 30;


%% Antenna Properties

Antenna = phased.IsotropicAntennaElement('FrequencyRange',[5e9 15e9]);

SensorMotion = phased.Platform('InitialPosition',[0; 0; 0], 'Velocity',[0; 0; 0]);

Radiator = phased.Radiator('Sensor',Antenna, 'OperatingFrequency', fc);

Collector = phased.Collector('Sensor',Antenna, 'OperatingFrequency',fc);

RCS_Reduction_Factor = 0.1; % takes value between 0 and 1

   

%% Targets
tgtpos = [[2024.66;0;0],[3518.63;0;0],[3845.04;0;0]]; %meters
tgtvel = [[100;100;0],[0;0;0],[0;0;0]];                % m/s
tgtmotion = phased.Platform('InitialPosition',tgtpos,'Velocity',tgtvel);

tgtrcs = [1.6 2.2 1.05];  % m²
target = phased.RadarTarget('MeanRCS',tgtrcs,'OperatingFrequency',fc);
%% Propagation Environment

channel = phased.FreeSpace('SampleRate',fs, 'TwoWayPropagation',true, 'OperatingFrequency',fc);

%% Signal Synthesis

fast_time_grid = unigrid(0,1/fs,1/PRF,'[)'); % Within each pulse (Colums)
slow_time_grid = (0:num_pulse_int-1)/PRF; % Time between pulses (Rows)

receiver.SeedSource = 'Property';
receiver.Seed = 2007;

% Pre-allocate array for improved processing speed
rxpulses = zeros(numel(fast_time_grid),num_pulse_int);

for m = 1:num_pulse_int

    % Update sensor and target positions
    [sensorpos,sensorvel] = SensorMotion(1/PRF);
    [tgtpos,tgtvel] = tgtmotion(1/PRF);

    % Calculate the target angles as seen by sthe sensor
    [tgtrng,tgtang] = rangeangle(tgtpos,sensorpos);

    % Simulate propagation of pulse in direction of targets
    pulse = waveform();
    [txsig,txstatus] = transmitter(pulse);
    txsig = Radiator(txsig,tgtang);
    txsig = channel(txsig,sensorpos,tgtpos,sensorvel,tgtvel);

    % Reflect pulse off of targets
    tgtsig = target(txsig);

    % Receive target returns at sensor
    rxsig = Collector(tgtsig,tgtang);
    rxpulses(:,m) = receiver(rxsig,~(txstatus>0));
end
%% Range-Doppler Response

RangeDopplerEx_MF_NFFTDOP = 1024;
RangeDopplerEx_MF_Fs = fs;
RangeDopplerEx_MF_Fc = fc;

response = phased.RangeDopplerResponse('DopplerFFTLengthSource','Property', ...
   'DopplerFFTLength',RangeDopplerEx_MF_NFFTDOP, ...
   'SampleRate',RangeDopplerEx_MF_Fs,'DopplerOutput','Speed', ...
   'OperatingFrequency',RangeDopplerEx_MF_Fc);

RangeDopplerEx_MF_X = rxpulses;
RangeDopplerEx_MF_Coeff = getMatchedFilter(waveform);
[resp,rng_grid,dop_grid] = response(RangeDopplerEx_MF_X, RangeDopplerEx_MF_Coeff);

figure,
imagesc(dop_grid,rng_grid,mag2db(abs(resp)));
xlabel('Speed (m/s)');
ylabel('Range (m)');
title('Range-Doppler Map Oerlikon Coherent');

%% --- Clutter Parametreleri ---
weibullShape = 2;      % Weibull shape (k)
weibullScale = 1;      % Weibull scale (lambda)
rho = 0.9;             % AR(1) korelasyon katsayısı

%% --- Döngü İçinde Clutter Oluştur ve Uygula ---
for m = 1:num_pulse_int
    % Sensor ve hedef konumlarını güncelle
    [sensorpos,sensorvel] = SensorMotion(1/PRF);
    [tgtpos,tgtvel] = tgtmotion(1/PRF);
    [tgtrng,tgtang] = rangeangle(tgtpos,sensorpos);

    % Radar darbesi oluştur
    pulse = waveform();
    [txsig,txstatus] = transmitter(pulse);
    txsig = Radiator(txsig, tgtang);
    txsig = channel(txsig, sensorpos, tgtpos, sensorvel, tgtvel);
    tgtsig = target(txsig);
    rxsig = Collector(tgtsig, tgtang);

    % --- Clutter üret (Weibull + AR + kompleks) ---
    N = numel(fast_time_grid);
    clutter = wblrnd(weibullScale, weibullShape, [N,1]);     % Weibull
    clutter = filter(1, [1 -rho], clutter);                  % AR(1)
    clutter = clutter .* exp(1j * 2 * pi * rand(N,1));       % Kompleks hale getir

    % --- Hedef sinyaline clutter ekle ve alıcıya gönder ---
    rxpulses(:,m) = receiver(rxsig + clutter, ~(txstatus > 0));
end
