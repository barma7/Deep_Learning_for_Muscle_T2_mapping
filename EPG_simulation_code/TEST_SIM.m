clear all; close all;

% Define folder for pulses
homeFldr = 'C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PHILIPS\PaperSharing\GeneratePulses_and_SliceProfiles\SINC_pulses\TWB2';
exc = readmatrix(fullfile(homeFldr, '90', 'SLR', 'pulse_profile_half.txt'),'delimiter',' ');
ref = readmatrix(fullfile(homeFldr, '180', 'SLR', 'pulse_profile_half.txt'),'delimiter',' ');

EchoSpacing = 7.6; %in ms (milliseconds)
ETL = 17;
% Set options for water
optw.esp = EchoSpacing;
optw.etl = ETL;
optw.mode = 's';
optw.RFe.alpha = exc'.*(pi/180);
optw.RFr.alpha = ref'.*(pi/180);
optw.T1 = 1400;
optw.Nz = size(ref,1);
optw.debug = 0;

% Set options for fat
optf.esp = EchoSpacing;
optf.etl = ETL;
optf.mode = 's';
optf.RFe.alpha = exc'.*(pi/180);
optf.RFr.alpha = ref'.*(pi/180);
optf.T1 = 400;
optf.Nz = size(ref,1);
optf.debug = 0;

%% Uniform Sampling
b1 = 0.8;  
T2_f = 150; % T2 fat in ms (milliseconds)
T2_w = 35; % T2 water in ms (milliseconds)
ff = 0.05;

tic
signal = cpmg_muscle_FSEsig(T2_w, T2_f, ff, b1, optw, optf);
toc

norm(signal);

tic
signal2 = cpmg_muscle_epg(ETL, EchoSpacing, exc, ref, b1, ff, T2_w, T2_f);
toc

taxis = EchoSpacing:EchoSpacing:ETL*EchoSpacing;
colors = {[55,126,184]./255, [228, 26, 28]./255, [77,175,74]./255, [255,127,0]./255,...
    [152,78,163]./255, [166,86,40]./255, [247,129,191]./255};

figure(); 
plot(taxis, abs(signal), "-o", "MarkerFaceColor",colors{1}, "LineWidth",1.5);
hold on;
plot(taxis, abs(signal2), "--", "MarkerFaceColor",colors{2}, "LineWidth",1.5)
hold on;
xlabel('time (ms)');
ylabel('Signal Intensity (a.u.)');
set(gca,'FontSize',14, 'Box','on','LineWidth', 2);

figure(); 
plot(taxis, abs(signal./norm(signal)), "-o", "MarkerFaceColor",colors{1}, "LineWidth",1.5);
hold on;
plot(taxis, abs(signal2./norm(signal2)), "--", "MarkerFaceColor",colors{2}, "LineWidth",1.5)
hold on;
xlabel('time (ms)');
ylabel('Signal Intensity (a.u.)');
set(gca,'FontSize',14, 'Box','on','LineWidth', 2);

N0 = 1e-3;
[y, SNR_g] = add_additive_noise_var((signal./norm(signal))',N0);

    