% CREATE SIMENS PULSE PROFILES FOR 2D MESE WITH GRADIENTS
clear all; close all;
% Define some specs
Nb_split = 100;
xm =  (Nb_split -1)/2;
dt90 = 2.0e-3; %[s]
dt180 = 2.0e-3; %[s]
tbw = 2;
thickness_factor = 1.2;
bw90 = tbw/dt90; % [Hz]
bw180 = tbw/dt180; %[Hz]

svFldr = fullfile(pwd, strcat("TWB",num2str(tbw)));
svFldr90 = fullfile(svFldr, "90");
svFldr180 = fullfile(svFldr, "180");
if ~exist(fullfile(svFldr90), 'dir') || ~exist(fullfile(svFldr180), 'dir')
    mkdir(fullfile(svFldr90))
    mkdir(fullfile(svFldr180))
end

dwell_90 = dt90/Nb_split;
dwell_180 = dt180/Nb_split;
time90 = dt90/Nb_split:dt90/Nb_split:dt90;
time180 = dt180/Nb_split:dt180/Nb_split:dt180;
x = (-xm:xm)./xm;
pulse_temp90 = round(hann(length(x))'.*sinc(x*tbw/2), 5);
pulse_temp180 = round(hann(length(x))'.*sinc(x*tbw/2), 5);

pulse90 = pulse_temp90;
pulse180 = pulse_temp180;

%figure; plot(time90, pulse90, time180, pulse180);

% Create Gradients according to slice thickness Gz = BW/(gamma*DZ)
gamma = 42.58*10^6;
z90 = 0.006; %[m]
z180 = z90*thickness_factor; %[m]
slwr = 180; %[T/m/s];

g90 = bw90/(gamma*z90);
g180 = bw180/(gamma*z180);

grad_90_slice_selective = g90.*ones(1,length(pulse_temp90));
grad_180_slice_selective = g180.*ones(1,length(pulse_temp180));

%figure; plot(time90, grad_90_slice_selective, time180, grad_180_slice_selective);

%% CREATE REPHASING LOBE AND CRUSCHING GRADIENTS FOR BLOCH SIMULATIONS
%Create rephasing gradient lobe
area = dt90/2 * g90;
Nb_pnts_reph = 30;
D_reph = dwell_90*Nb_pnts_reph;
greph = area/(dwell_90*Nb_pnts_reph);
grad_90_rephasing = -greph*ones(1, Nb_pnts_reph);

% Evaluate rephasing 90 and Crushing Gradient 180
Acrusher = (4*pi)/(2*pi*gamma*z180);
Nb_pnts_crush = round(Nb_pnts_reph/2);
D_crusher = dwell_180*Nb_pnts_crush;
gmax = Acrusher/D_crusher;

grad_180_curhser = gmax.*ones(1, Nb_pnts_crush);

% Putting things togheter
grad_90_total = cat(2, grad_90_slice_selective, grad_90_rephasing);
grad_180_total = cat(2, grad_180_curhser, grad_180_slice_selective, grad_180_curhser);

pulse90_total = cat(2, pulse90, zeros(1, length(grad_90_rephasing)));
pulse180_total = cat(2, zeros(1, length(grad_180_curhser)), pulse180, zeros(1, length(grad_180_curhser)));

time90_temp = dwell_90:dwell_90:dwell_90*length(pulse90_total);
time180_temp = dwell_180:dwell_180:dwell_180*length(pulse180_total);

figure; plot(time90_temp, grad_90_total, time180_temp, grad_180_total);
figure; plot(time90_temp, pulse90_total, time180_temp, pulse180_total);
% Spec table
colnames = {'FA', 'Pulse Lenght (ms)', 'Pulse_Nb_points', 'TBW', 'Slice Thickness (mm)',...
    'thickness factor', 'dwell (s)', 'SS g (T/m)', 'Reph/crush g (T/m)',...
    'Reph/crush time (ms)'};
T = table([90; 180], [dt90; dt180].*1000, [Nb_split; Nb_split], [tbw; tbw],...
    [z90; z180].*1000, [NaN; thickness_factor],...
    [dwell_90; dwell_180], [g90; g180], [greph; gmax],...
    [D_reph; D_crusher].*1000, 'VariableNames', colnames);

% save waveform for 90
writematrix(pulse90_total, fullfile(svFldr90,'Pulse_total.txt'),'Delimiter',',');
writematrix(grad_90_total, fullfile(svFldr90,'Grad_total.txt'), 'Delimiter',',');
writematrix(pulse90, fullfile(svFldr90,'Pulse_only.txt'),'Delimiter',',');
writematrix(grad_90_slice_selective, fullfile(svFldr90,'Grad_only.txt'), 'Delimiter',',');
writematrix(D_reph*1000, fullfile(svFldr90, 'DeltaTime_rephaser_ms.txt'));
writematrix(dwell_90, fullfile(svFldr90, 'Dwell_pulse_s.txt'));
writematrix(dt90*1000, fullfile(svFldr90, 'DeltaTime_ms.txt'));
% save waveform for 180
writematrix(pulse180_total, fullfile(svFldr180,'Pulse_total.txt'),'Delimiter',',');
writematrix(grad_180_total, fullfile(svFldr180, 'Grad_total.txt'), 'Delimiter',',');
writematrix(pulse180, fullfile(svFldr180,'Pulse_only.txt'),'Delimiter',',');
writematrix(grad_180_slice_selective, fullfile(svFldr180,'Grad_only.txt'), 'Delimiter',',');
writematrix(dt180*1000, fullfile(svFldr180, 'DeltaTime_ms.txt'));
writematrix(D_crusher*1000, fullfile(svFldr180, 'DeltaTime_crusher_ms.txt'));
writematrix(dwell_180, fullfile(svFldr180, 'Dwell_pulse_s.txt'));

% save specs
writetable(T, fullfile(svFldr, 'PulseSpecs.csv'));