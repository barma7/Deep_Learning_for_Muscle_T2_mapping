%% CREATE CPMG Examples
clear all; close all;

% IMPORT PULSE PROFILES
homeFldr = 'C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PHILIPS\PaperSharing\GeneratePulses_and_SliceProfiles\SINC_pulses\TWB2';
exc = readmatrix(fullfile(homeFldr, '90', 'SLR', 'pulse_profile_half.txt'),'delimiter',' ');
ref = readmatrix(fullfile(homeFldr, '180', 'SLR', 'pulse_profile_half.txt'),'delimiter',' ');

% Define Sving fulders
SvFldr_home = 'C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PHILIPS\PaperSharing\sim-data';
exp_fldr = 'cpmg_ET_7.6ms_ETL_17\100000-examples_t2w-10-110-t2f-50-200-b1-1.2';

if ~exist(fullfile(SvFldr_home, exp_fldr), 'dir')
    mkdir(fullfile(SvFldr_home, exp_fldr))
end
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
Nb_exmpls = 100000;
rng('shuffle');
lowT2_w = 10;
highT2_w = 110;
lowT2_f = 50;
highT2_f = 200;
B1min = 0.3;
B1max = 1.2;

rng(5); % to assure reproducibility (optional)

T2w = lowT2_w + (highT2_w-lowT2_w).*rand(Nb_exmpls,1);
T2f = lowT2_f + (highT2_f-lowT2_f).*rand(Nb_exmpls,1);
B1s = B1min + (B1max-B1min).*rand(Nb_exmpls,1);
LUT = cat(2,T2w, T2f, B1s);

dictionary_water = zeros(Nb_exmpls,ETL);
dictionary_fat = zeros(Nb_exmpls,ETL);

tic
parfor i = 1:Nb_exmpls
    [dictionary_water(i,:), dictionary_fat(i, :)] =  cpmg_muscle_FSEsig2(T2w(i), T2f(i), B1s(i), optw, optf);
end
toc

writematrix(LUT, fullfile(SvFldr_home,exp_fldr,'LUT.txt'));
writematrix(ETL, fullfile(SvFldr_home,exp_fldr,'ETL.txt'));
writematrix(exc, fullfile(SvFldr_home,exp_fldr,'exc_profile.txt'));
writematrix(ref, fullfile(SvFldr_home,exp_fldr,'ref_profile.txt'));

save(fullfile(SvFldr_home,exp_fldr,'dictionary_water.mat'),'dictionary_water')
save(fullfile(SvFldr_home,exp_fldr,'dictionary_fat.mat'),'dictionary_fat')