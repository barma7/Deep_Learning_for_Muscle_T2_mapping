close all; clear all;
%RF tool path (NB: Change path accordingly)
rftool_path = "C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Matlab_miscellaneous\Spectral-Spatial-RF-Pulse-Design\rf_tools";
addpath(rftool_path, fullfile(rftool_path, "mex_files"));

home_fldr = fullfile(convertCharsToStrings(pwd), "TWB2");
fldrs = ["90", "180"];

Npnts = 128;

for i = 1:length(fldrs)
    pulse_folder = fullfile(home_fldr, fldrs(i));
    pulse = readmatrix(fullfile(pulse_folder, "Pulse_only.txt"), 'Delimiter',',');
    fa = str2double(fldrs(i)).*pi/180;
    pulse = fa.*(pulse./sum(pulse));
    dwell = readmatrix(fullfile(pulse_folder, "Dwell_pulse_s.txt")).*1000;
    grad = readmatrix(fullfile(pulse_folder, "Grad_only.txt"));

    svFldr = fullfile(pulse_folder,"SLR");
    if ~exist(fullfile(svFldr), 'dir')
        mkdir(fullfile(svFldr))
    end

    tpulse = dwell*length(pulse);
    Dx = 10;
    x = linspace(-Dx/2, Dx/2,Npnts);
    xint = linspace(-.5, .5, Npnts);

    if i == 1
        g = grad(round(length(pulse)/2))*100;
        ab = abr(pulse,x);
        sp = abs(ab2ex(ab))./max(abs(ab2ex(ab)));
        xs = gt2cm(x,g,tpulse);
        spint = interp1(xs, sp, xint);
        spint(isnan(spint)) = 0;
    else
        g = grad(round(length(pulse)/2))*100;
        ab = abr(pulse,x);
        sp = abs(ab2se(ab))./max(abs(ab2se(ab)));
        xs = gt2cm(x,g,tpulse);
        spint = interp1(xs, sp, xint);
        spint(isnan(spint)) = 0;
    end
    
    writematrix(str2double(fldrs(i)).*spint', fullfile(svFldr, "pulse_profile.txt"));
    writematrix(str2double(fldrs(i)).*spint(Npnts/2:end)', fullfile(svFldr, "pulse_profile_half.txt"));
    figure(1);
    plot(xint, fa.*abs(spint));
    hold on;
    figure(2);
    plot(fa.*abs(spint));
    hold on;
end


rmpath(rftool_path, fullfile(rftool_path, "mex_files"));