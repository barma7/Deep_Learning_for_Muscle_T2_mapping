clear all; close all;

home_path = "C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PaperSharing\DATA\SIEMENS";

subjects = ["01", "02", "03", "04", "05"];

ffmap_names = ["ff_s1.nii"];
t2map_names = ["t2_s1.nii"];
b1map_names = ["b1_s1.nii"];
rsq_names = ["rsq_s1.nii"];
snr_names = ["snr_s1.nii"];
tse_names = ["mese_vol.nii.gz"];

sv_name = "MAPS/dictionary/SINC/TBW2";

fat_calibration = 1;
fat_cal_name = "mese_vol.nii.gz";

%% PATHS FOR BLOCH SIMULATION FORMATION
% Add EPG simulation to MATLAB path
epg_path = "C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PaperSharing\EPG_simulation_codes";
addpath(epg_path);

% load slice profiles
homeFldr = 'C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PaperSharing\GeneratePulses_and_SliceProfiles\SINC_pulses\TWB2';
exc = readmatrix(fullfile(homeFldr, '90', 'SLR', 'pulse_profile_half.txt'),'delimiter',' ');
ref = readmatrix(fullfile(homeFldr, '180', 'SLR', 'pulse_profile_half.txt'),'delimiter',' ');

ETL = 17;
TE = 7.5;

% Set options for water
optw.esp = TE;
optw.etl = ETL;
optw.mode = 's';
optw.RFe.alpha = exc'.*(pi/180);
optw.RFr.alpha = ref'.*(pi/180);
optw.T1 = 1400;
optw.Nz = size(ref,1);
optw.debug = 0;

% Set options for fat
optf.esp = TE;
optf.etl = ETL;
optf.mode = 's';
optf.RFe.alpha = exc'.*(pi/180);
optf.RFr.alpha = ref'.*(pi/180);
optf.T1 = 400;
optf.Nz = size(ref,1);
optf.debug = 0;

if fat_calibration == 1
    disp("Creating dictionary for fat calibration...")
    % CREATE DICTIONARY FOR FAT CALIBRATION 
    T2ws = 10:0.5:80;
    T2fs = 50:0.25:250;
    B1s = 0.5:0.025:1.2;
    [dictionary_fat, LUT_fatCal] = create_fat_dictionary(optw, optf, T2ws, T2fs, B1s, 0.9);
end
   
for k = 1:length(subjects)
    disp(strcat("Processing Subject ", subjects(k)))
    
    sub = subjects(k);
    SvFldr = fullfile(home_path, sub, sv_name);

    if ~exist(SvFldr, 'dir')
        mkdir(SvFldr);
    end
    
    data_fldr = fullfile(home_path, sub);
    
    %check if FatCal is present in folder
    if fat_calibration == 0
        
        T = readmatrix(fullfile(SvFldr, "fat_calibration.csv"));
        t2fat = T(2,3);
        
    else
        % calibrate T2fat
        disp("pefroming fat Calibration...")
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
        
        % Remove small connected components
        for sl = round(nslice / 2) - 2:round(nslice / 2) + 2
            fat_mask(:,:,sl) = bwareaopen(fat_mask(:,:,sl), 40);
        end

        data_fat_flat = reshape(Data_fat, nrow*ncol*nslice, necho);
        fat_mask_flat = fat_mask(:);
        data_fat = data_fat_flat(fat_mask_flat > 0, :);
        
        %normalize data
        for i = 1:size(data_fat,1)
            data_fat(i,:) = data_fat(i,:)/norm(data_fat(i,:));
        end
       
        %normalize dictioanry
        for i = 1:size(dictionary_fat,1)
            dictionary_fat(i,:) = abs(dictionary_fat(i,:))/norm(abs(dictionary_fat(i,:)));
        end
        
        [T2fat, T2water, B1fit, tElapsed] = FatCal_dictionary_matching_dot(LUT_fatCal, dictionary_fat, data_fat);
        disp(['Fitted Fat data in : ', num2str(tElapsed)]);

        % Take only pixels where b1 is grater than 0.7
        valid_pixels = find(B1fit > 0.7);
        T2water = T2water(valid_pixels);
        T2fat = T2fat(valid_pixels);
        B1fit = B1fit(valid_pixels);
        
        % plot 2D histogram
        figure;
        hist3(cat(2, T2water, T2fat),[10, 10], 'CdataMode','auto')
        xlabel('T2 water (ms)')
        ylabel('T2 fat (ms)')
        colorbar
        view(2)
       
        t2_fat_mean = mean(T2fat, 'omitnan');
        t2_water_mean = mean(T2water, 'omitnan');

        t2_fat_geom = geomean(T2fat, 'omitnan');
        t2_water_geom = geomean(T2water, 'omitnan');
        
        % Save Data
        writematrix(T2water, fullfile(SvFldr, "SCF_t2_water.csv"));
        writematrix(T2fat, fullfile(SvFldr, "SCF_t2_fat.csv"));
        
        Tissue = ["Water";"Fat"];
        t2mean = cat(1, t2_water_mean, t2_fat_mean);
        t2geom = cat(1, t2_water_geom, t2_fat_geom);

        t2fat = t2_fat_mean;

        T = table(Tissue, t2mean, t2geom);
        writetable(T, fullfile(SvFldr, "fat_calibration.csv"));

        disp("Fat Calibration Completed. T2 fat is");
        disp(t2fat);
    end
    disp("Building Dictionary for T2 mapping...")
    
    T2ws = 10:0.2:80;
    T2fs = t2fat;
    B1s = 0.5:0.02:1.2;
    FFs = 0.01:0.015:1;
    
    [dictionary, LUT] = create_muscle_dictionary(optw, optf, T2ws, T2fs, B1s, FFs);

    % Normalize Dictionary
    for i=1:length(dictionary(:,1))
        dictionary(i,:) = dictionary(i,:)./norm(dictionary(i,:));
    end

    tic 
    for stack = 1:length(tse_names)
        % Create FF, T2 and B1 maps

        data_nifti = single(niftiread(fullfile(data_fldr,tse_names(stack))));
        info_nifti = niftiinfo(fullfile(data_fldr,tse_names(stack)));
        
        nb_row = size(data_nifti,1);
        nb_col = size(data_nifti,2);
        nb_slices = size(data_nifti,3);
        ETL = size(data_nifti,4);
        
        FFmap = single(zeros(nb_row, nb_col, nb_slices));
        T2map = single(zeros(nb_row, nb_col, nb_slices));
        B1map = single(zeros(nb_row, nb_col, nb_slices));
        
        tstart = tic;
        for slc = 1:nb_slices
            
            FF = zeros(nb_row * nb_col,1);
            T2 = zeros(nb_row * nb_col,1);
            B1 = zeros(nb_row * nb_col,1);
            
            data_temp = squeeze(data_nifti(:, :, slc, :));
            data_temp = reshape(data_temp, nb_row*nb_col, ETL);

            non_zero_indexes = find(data_temp(:, 1) > 50);
            data_temp = data_temp(non_zero_indexes, :);
            
            % Normalize data
            for i=1:length(data_temp(:,1))
                data_temp(i,:) = data_temp(i,:)./norm(data_temp(i,:));
            end
            
            [FFfit, T2fit, B1fit, tElapsed1] = dictionary_matching_dot(LUT(:,[1,3,4]), dictionary, data_temp);
            
            FF(non_zero_indexes,1) = FFfit;
            T2(non_zero_indexes,1) = T2fit;
            B1(non_zero_indexes,1) = B1fit;

            FFmap(:, :, slc) = reshape(FF, nb_row, nb_col);
            T2map(:, :, slc) = reshape(T2, nb_row, nb_col);
            B1map(:, :, slc) = reshape(B1, nb_row, nb_col);
        end
        tElapsed = toc(tstart);
        disp(['Fitted non zero element in subject in : ', num2str(tElapsed)]);

        % save maps
        info_out = info_nifti;
        info_out.ImageSize = info_nifti.ImageSize(1:3);
        info_out.Datatype = 'single';
        info_out.PixelDimensions = info_nifti.PixelDimensions(1:3);
        
        niftiwrite(FFmap, fullfile(SvFldr, ffmap_names(stack)), info_out);
        niftiwrite(T2map, fullfile(SvFldr, t2map_names(stack)), info_out);
        niftiwrite(B1map, fullfile(SvFldr, b1map_names(stack)), info_out);
        
    end
    toc
end
    
        
        
        




