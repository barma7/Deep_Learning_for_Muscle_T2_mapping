%% Script to perform estimation
clear all; close all;

subjects = ["AM1705", "AM1708", "AM1709", "AM1714", "AM1716"];

ffmap_names = ["ff_s1.nii", "ff_s2.nii", "ff_s3.nii"];
t2map_names = ["t2_s1.nii", "t2_s2.nii", "t2_s3.nii"];
b1map_names = ["b1_s1.nii", "b1_s3.nii", "b1_s3.nii"];
rsq_names = ["rsq_s1.nii", "rsq_s2.nii", "rsq_s3.nii"];
snr_names = ["snr_s1.nii", "snr_s2.nii", "snr_s3.nii"];
tse_names = ["TSE_s1.nii.gz", "TSE_s2.nii.gz", "TSE_s3.nii.gz"];

sv_name = "MAPS/SNLSQ/SINC/TBW2";

fat_calibration = 1;
fat_cal_name = "TSE_s2.nii.gz";

%% PATHS FOR BLOCH SIMULATION FORMATION
% Add EPG simulation to MATLAB path
epg_path = "C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PaperSharing\EPG_simulation_codes";
addpath(epg_path);

% load slice profiles
homeFldr = 'C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PaperSharing\GeneratePulses_and_SliceProfiles\SINC_pulses\TWB2';
exc = readmatrix(fullfile(homeFldr, '90', 'SLR', 'pulse_profile_half.txt'),'delimiter',' ');
ref = readmatrix(fullfile(homeFldr, '180', 'SLR', 'pulse_profile_half.txt'),'delimiter',' ');

ETL = 17;
TE = 7.6; % [ms]

% Set options for water
opt.optw.esp = TE;
opt.optw.etl = ETL;
opt.optw.mode = 's';
opt.optw.RFe.alpha = exc'.*(pi/180);
opt.optw.RFr.alpha = ref'.*(pi/180);
opt.optw.T1 = 1400;
opt.optw.Nz = size(ref,1);
opt.optw.debug = 0;

% Set options for fat
opt.optf.esp = TE;
opt.optf.etl = ETL;
opt.optf.mode = 's';
opt.optf.RFe.alpha = exc'.*(pi/180);
opt.optf.RFr.alpha = ref'.*(pi/180);
opt.optf.T1 = 400;
opt.optf.Nz = size(ref,1);
opt.optf.debug = 0;

%   Numeric fitting options
opt.lsq.fopt = optimset('lsqcurvefit');
opt.lsq.fopt.TolX = 1e-3;     %   Fitting accuracy: 0.1 ms
opt.lsq.fopt.TolFun = 1.0e-9;
opt.lsq.fopt.MaxIter = 400;
opt.lsq.fopt.Display = 'off';

list_sub_id = zeros(length(subjects),1);
list_time_fat_cal = zeros(length(subjects),1);
list_time_water = zeros(length(subjects),1);

%% FIT of DATA
for k = 1:length(subjects)
    
    list_sub_id(k) = k;
    sub = subjects(k);
    SvFldr = fullfile(home_path, sub, sv_name);

    if ~exist(SvFldr, 'dir')
        mkdir(SvFldr);
    end
    
    data_fldr = fullfile(home_path, sub);
    
    if fat_calibration == 1
        %% FAT CAL FITTING OPTIONS
        % This is T2 fat and T2 water separation
        opt.lsq.Icomp.X0   = [30, 155, 1];      %   Starting point (1 x 3) [T2w, T2f, B1]
        opt.lsq.Icomp.XU   = [80, 250, 1.2];  %   Upper bound (1 x 3)
        opt.lsq.Icomp.XL   = [10, 50, 0.4];  %   Lower bound (1 x 3)
        
        % LOAD DATA
        Data_fat = niftiread(fullfile(data_fldr, fat_cal_name));

        nrow = size(Data_fat, 1);
        ncol = size(Data_fat, 2);
        nslice = size(Data_fat, 3);
        necho = size(Data_fat, 4);

        % Cluster Segmentation & select data from subcutaneous fat
        fat_mask = zeros(nrow, ncol, nslice);
        temp = single(Data_fat(:, :, round(nslice / 2) - 2:round(nslice / 2) + 2, end));
        clustered = imsegkmeans3(temp,2);
        clustered_flat = clustered(:);
        s1 = length(clustered_flat(clustered_flat == 1));
        s2 = length(clustered_flat(clustered_flat == 2));
        
        if s1 >= s2
            clustered(clustered < 2) = 0;
            clustered(clustered > 0) = 1;
        else
            clustered(clustered > 1) = 0;
        end  

        fat_mask(:, :, round(nslice / 2) - 2:round(nslice / 2) + 2) = clustered; 

        data_fat_flat = reshape(Data_fat, nrow*ncol*nslice, necho);
        fat_mask_flat = fat_mask(:);
        data_fat = data_fat_flat(fat_mask_flat > 0, :);

        % PERFORM FITTING
        Nb_examples = size(data_fat,1);

        x_est = zeros(Nb_examples, 3);
    
        ub = opt.lsq.Icomp.XU;
        lb = opt.lsq.Icomp.XL;

        x0 = opt.lsq.Icomp.X0;

        rhos = zeros(Nb_examples,1);

        tstart = tic;
        parfor i=1:Nb_examples

            ydata = data_fat(i,:)';
        
            %figure(2)
            %plot(ydata,"LineWidth",2,"Color","black"); hold on;
            
            % Inizialize starting point
            F = @(x, varargin)fit_FAT_EPG(x, ydata, opt);
        
            x_est(i,:)= lsqcurvefit(F,x0, opt, ydata, lb, ub, opt.lsq.fopt);
            
            %Estimate rhos
            rhos(i) = fit_linear_coefficient_FAT_EPG(x_est(i,:), ydata, opt);
            %close(2);
       
        end
        tElapsed = toc(tstart)

        list_time_fat_cal(k) = tElapsed;
        
        t2_water_scatter = x_est(:,1);
        t2_fat_scatter = x_est(:,2);
        m0fat = rhos(:);

        figure;
        hist3(cat(2, t2_water_scatter, t2_fat_scatter),[10, 10], 'CdataMode','auto')
        xlabel('T2 water (ms)')
        ylabel('T2 fat (ms)')
        colorbar
        view(2)
       
        t2_fat_mean = mean(t2_fat_scatter, 'omitnan');
        t2_water_mean = mean(t2_water_scatter, 'omitnan');

        t2_fat_geom = geomean(t2_fat_scatter, 'omitnan');
        t2_water_geom = geomean(t2_water_scatter, 'omitnan');
        
        % Save Data

        writematrix(t2_fat_scatter, fullfile(SvFldr, "SCF_t2_fat.csv"));
        writematrix(t2_water_scatter, fullfile(SvFldr, "SCF_t2_water.csv"));

        Tissue = ["Water";"Fat"];
        t2mean = cat(1, t2_water_mean, t2_fat_mean);
        t2geom = cat(1, t2_water_geom, t2_fat_geom);

        T2fat = t2_fat_geom;

        T = table(Tissue, t2mean, t2geom);
        writetable(T, fullfile(SvFldr, "fat_calibration.csv"));

    else
        T = readmatrix(fullfile(SvFldr, "fat_calibration.csv"));
        T2fat = T(2,3);
    end

    %% FIT Muscles

    % Fixed T2fat, fit FF, T2water and B1
    opt.lsq.Icomp.X0   = [30, 1];      %   Starting point (1 x 2) [T2w, B1]
    opt.lsq.Icomp.XU   = [80, 1.2];  %   Upper bound (1 x 2)
    opt.lsq.Icomp.XL   = [10, 0.4];  %   Lower bound (1 x 2)
    
    opt.T2fat = T2fat;
    
    for stack = 1:length(tse_names)

        % LOAD DATA
        Data = niftiread(fullfile(data_fldr, tse_names(stack)));
        info = niftiinfo(fullfile(data_fldr, tse_names(stack)));

        nrow = size(Data, 1);
        ncol = size(Data, 2);
        nslice = size(Data, 3);
        necho = size(Data, 4);

        FFmap = zeros(nrow*ncol*nslice, 1);
        T2map = zeros(nrow*ncol*nslice, 1);
        B1map = zeros(nrow*ncol*nslice, 1);
        R2map = zeros(nrow*ncol*nslice, 1);

        % Flat Data and select MAG > 50
        Data_flat = reshape(Data, nrow*ncol*nslice, necho);
        non_zero_indexes = find(Data_flat(:,2) > 50);
        data = Data_flat(non_zero_indexes, :);
        
        % PERFORM FITTING
        Nb_examples = size(data,1);

        x_est = zeros(Nb_examples, 2);
    
        ub = opt.lsq.Icomp.XU;
        lb = opt.lsq.Icomp.XL;

        x0 = opt.lsq.Icomp.X0;

        rhos = zeros(Nb_examples,2);
        Rsq = zeros(Nb_examples,1);

        tstart = tic;
        parfor i=1:Nb_examples
            ydata = data(i,:)';

            %figure(2)
            %plot(ydata,"LineWidth",2,"Color","black"); hold on;
            
            % Inizialize starting point
            F = @(x, varargin)fit_MUSCLE_EPG(x, ydata, opt);
        
            x_est(i,:)= lsqcurvefit(F,x0, opt, ydata, lb, ub, opt.lsq.fopt);
            
            %Estimate rhos
            [rhos(i,:), Rsq(i)] = fit_linear_coefficient_MUSCLE_EPG(x_est(i,:), ydata, opt);
            %close(2);
       
        end
        tElapsed = toc(tstart)
        list_time_water(k) = tElapsed;
        
        T2map(non_zero_indexes) = x_est(:, 1);
        B1map(non_zero_indexes) = x_est(:, 2);
        FFmap(non_zero_indexes) = 100.*(abs(rhos(:,2))./sum(abs(rhos),2));
        R2map(non_zero_indexes) = Rsq;
        
        T2map = reshape(T2map, nrow, ncol, nslice);
        B1map = reshape(B1map, nrow, ncol, nslice);
        FFmap = reshape(FFmap, nrow, ncol, nslice);
        R2map = reshape(R2map, nrow, ncol, nslice);

        % save maps
        info_out = info;
        info_out.ImageSize = info.ImageSize(1:3);
        info_out.Datatype = 'double';
        info_out.PixelDimensions = info.PixelDimensions(1:3);

        niftiwrite(T2map, fullfile(SvFldr, t2map_names(stack)),info_out);
        niftiwrite(FFmap, fullfile(SvFldr, ffmap_names(stack)),info_out);
        niftiwrite(B1map, fullfile(SvFldr, b1map_names(stack)),info_out);
        niftiwrite(R2map, fullfile(SvFldr, rsq_names(stack)),info_out);
    end
end

colnames = {'sub_id', 'time fat cal', 'time water fit'};
T = table(list_sub_id, list_time_fat_cal, list_time_water, 'VariableNames',colnames);

writetable(T, fullfile(home_path, 'processing_log_EPG_NLSQ_PHILIPS.csv'));








