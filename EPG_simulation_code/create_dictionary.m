%% CREATE CPMG Examples
clear all; close all;

% IMPORT PULSE PROFILES
homeFldr = 'C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PHILIPS\PaperSharing\GeneratePulses_and_SliceProfiles\SINC_pulses\TWB2';
exc = readmatrix(fullfile(homeFldr, '90', 'SLR', 'pulse_profile_half.txt'),'delimiter',' ');
ref = readmatrix(fullfile(homeFldr, '180', 'SLR', 'pulse_profile_half.txt'),'delimiter',' ');

SvFldr_home = 'C:\ShareFldr-Ubuntu\Lab-work\Code\Projects\Muscle-WaterT2\SIMENS\sim-data';
exp_fldr = 'cpmg_ET_7.5ms_ETL_17\dictionary_t2w-80-t2f-155-b1-1.2';

if ~exist(fullfile(SvFldr_home, exp_fldr), 'dir')
    mkdir(fullfile(SvFldr_home, exp_fldr))
end

EchoSpacing = 7.5; %in ms (milliseconds)
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
T2ws = 10:0.2:80;
T2fs = 155;
B1s = 0.5:0.01:1.2;
    
Nb_exmpls = length(T2ws)*length(T2fs)*length(B1s);

LUT = zeros(Nb_exmpls,2);
count = 1;
for k=1:length(T2ws)
    for l = 1: length(B1s)
        LUT(count,:) = cat(2,T2ws(k), B1s(l));
        count = count+1;
    end
end

dictionary_water = zeros(Nb_exmpls,ETL);
dictionary_fat = zeros(Nb_exmpls,ETL);

T2ws = LUT(:,1);
B1s = LUT(:,2);

tic
parfor i = 1:Nb_exmpls
    [dictionary_water(i,:), dictionary_fat(i, :)] =  cpmg_muscle_FSEsig2(T2ws(i), T2fs, B1s(i), optw, optf);
end
toc

writematrix(LUT, fullfile(SvFldr_home,exp_fldr,'LUT.txt'));
writematrix(ETL, fullfile(SvFldr_home,exp_fldr,'ETL.txt'));
writematrix(exc, fullfile(SvFldr_home,exp_fldr,'exc_profile.txt'));
writematrix(ref, fullfile(SvFldr_home,exp_fldr,'ref_profile.txt'));

%% Combine dictionary for FFs
FFs = 0:0.01:1;

dictionary = zeros(Nb_exmpls*length(FFs),ETL);
LUT_combined = zeros(Nb_exmpls*length(FFs), 3);

cnt = 1;
for i = 1:length(FFs)
    for j = 1:Nb_exmpls
        dictionary(cnt,:) = (1 - FFs(i)).*dictionary_water(j,:) + FFs(i).*dictionary_fat(j,:);
        LUT_combined(cnt,:) = cat(2, FFs(i), LUT(j,:));
        cnt = cnt + 1;
    end
end
       

save(fullfile(SvFldr_home,exp_fldr,'LUT_combined.mat'),'LUT_combined');
save(fullfile(SvFldr_home,exp_fldr,'dictionary.mat'),'dictionary');
    