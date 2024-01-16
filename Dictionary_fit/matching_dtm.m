%% matching experimental fingerprints with dictionary
clear all; close all;
Nb_pulses = 17;
apply_mask = 0;
phase_aligment = 0;
add_noise = 0;



% LOAD DICTIONARY 
dictFolder = '/bmrNAS/people/barma7/MATLAB_drive_vigata/Lab-work/2021/FWF_muscle_project/sim-data/marathon_study-cpmg_ET_7.6ms/correct_pulses/dictionary_b1_1.4/';
load(strcat(dictFolder,'dictionary.mat'));
dim = size(dictionary);
Dict = permute(dictionary,[1 3 2]);
Dict = reshape(Dict,[],size(dictionary,2),1);
load(strcat(dictFolder,'LUT_combined.mat'));
Dict_NMRpar_ = permute(LUT_combined,[1 3 2]);
Dict_NMRpar_ = reshape(Dict_NMRpar_,[],size(LUT_combined,2),1);
Dict_t1fat = LUT_combined(:,2:end,1);

n = 155;
   
[val,idx]=min(abs(Dict_t1fat(:,1)-n));
minVal=max(Dict_t1fat(idx,1));
Dict_NMRpar = cat(2, Dict_NMRpar_(Dict_NMRpar_(:,2) == minVal,1), Dict_NMRpar_(Dict_NMRpar_(:,2) == minVal,3:4));
sub_dictionary = Dict(Dict_NMRpar_(:,2) == minVal, :);

%######################## LOAD DATA SET ###################################
home_data_foldr = '/bmrNAS/people/barma7/MATLAB_drive_vigata/Lab-work/2021/FWF_muscle_project/';
dataFolder = strcat(home_data_foldr,'marathon_dataset/AM1709/');
data_dicom = dicomread(strcat(dataFolder,'IM_0013'));
save_folder = strcat(dataFolder,'dictionary_marco_correct_pulses/S3/b1_1.4/');
if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end

nb_row = size(data_dicom,1);
nb_col = size(data_dicom,2);

ETL = 17;
TE = 7.6;
tAxis = TE:TE:TE*ETL;
Nb_slices = round(size(data_dicom,4)/ETL);
Data = zeros(size(data_dicom,1), size(data_dicom,2), Nb_slices, ETL);
indx = 1;
for slice = 1:Nb_slices
    for echo = 1:ETL
        if (slice == Nb_slices && echo == ETL)
            Data(:,:,slice,echo) = squeeze(data_dicom(:,:,1,indx-1));
        else
            Data(:,:,slice,echo) = squeeze(data_dicom(:,:,1,indx));
            indx = indx + 1;
        end
    end
end
    
clear data dictionary;

%% PART II: USING THE COSINE DISTANCE TO PERFORM MATCHING
Nb_tests = nb_row*nb_col;

Dict_nor = single(zeros(size(sub_dictionary)));
for i=1:length(Dict_nor(:,1))
    Dict_nor(i,:) = sub_dictionary(i,:)./norm(sub_dictionary(i,:));
end


clear Dict

Divisore = 1000;
% break the Nb_tests in block of 10000 test data
if (Nb_tests/Divisore - round(Nb_tests/Divisore)) < 0.5 && (Nb_tests/Divisore - round(Nb_tests/Divisore))>0
    Nb_blocks = round(Nb_tests/Divisore)+1;
else
    Nb_blocks = round(Nb_tests/Divisore);
end


FFmap = zeros(nb_row, nb_col,Nb_slices);
T2map = zeros(nb_row, nb_col,Nb_slices);
B1map = zeros(nb_row, nb_col,Nb_slices);
Sim_data = complex(zeros(nb_row, nb_col,Nb_slices, ETL));

% Prepare MAP 

tstart = tic;
for s = 1:Nb_slices
    s
    %NORMALIZE DATA
    Data_temp = reshape(squeeze(Data(:,:,s,:)), nb_row*nb_col,ETL);
    Data_n = single(zeros(size(Data_temp)));
    for i=1:length(Data_n(:,1))
        Data_n(i,:) = Data_temp(i,:)./norm(Data_temp(i,:));
    end

    indx = 1;
    stop_indx = Divisore;
    start_indx = 1;
    FF = zeros(Nb_tests,1);
    T2 = zeros(Nb_tests,1);
    B1 = zeros(Nb_tests,1);
    Sim_sig = complex(zeros(Nb_tests,ETL));
    
    tic
    for b = 1:Nb_blocks
        if (b == Nb_blocks && stop_indx>Nb_tests)
            stop_indx = Nb_tests-(stop_indx-Divisore);
            dotMatrix = Data_n(start_indx:start_indx+stop_indx-1,:)*Dict_nor';    
            [maxDot_array, Mtchd_indx_array] = max(abs(dotMatrix),[],2);

            for j=1:length(maxDot_array) 
                Mtch_indx = Mtchd_indx_array(j);
                FF(indx) = Dict_NMRpar(Mtch_indx,1);
                T2(indx) = Dict_NMRpar(Mtch_indx,2);
                B1(indx) = Dict_NMRpar(Mtch_indx,3);
                Sim_sig(indx, :) = Dict_nor(Mtch_indx, :);

                indx = indx+1;
            end
        else
            dotMatrix = Data_n(start_indx:stop_indx,:)*Dict_nor';    
            [maxDot_array, Mtchd_indx_array] = max(abs(dotMatrix),[],2);

            for j=1:length(maxDot_array) 
                Mtch_indx = Mtchd_indx_array(j);
                FF(indx) = Dict_NMRpar(Mtch_indx,1);
                T2(indx) = Dict_NMRpar(Mtch_indx,2);
                B1(indx) = Dict_NMRpar(Mtch_indx,3);
                Sim_sig(indx, :) = Dict_nor(Mtch_indx, :);
                
                indx = indx+1;
            end
            start_indx = start_indx+Divisore;
            stop_indx = stop_indx+Divisore;
        end
    end
    toc
    FFmap(:,:,s) = reshape(FF, nb_row, nb_col);
    T2map(:,:,s) = reshape(T2, nb_row, nb_col);
    B1map(:,:,s) = reshape(B1, nb_row, nb_col);
    Sim_data(:,:,s,:) = reshape(Sim_sig, nb_row, nb_col, ETL);
    
end
tElapsed = toc(tstart)

%% Restore 3D maps
niftiwrite(FFmap,strcat(save_folder,'ff.nii'));
niftiwrite(T2map,strcat(save_folder,'t2.nii'));
niftiwrite(B1map,strcat(save_folder,'b1.nii'));
niftiwrite(Sim_data,strcat(save_folder,'Syntetic_data.nii'));

sl = 10;
figure;
imagesc(FFmap(:,:,sl)); caxis([0 1]); colorbar; colormap('hot');

figure;
imagesc(T2map(:,:,sl)); caxis([25 35]); colorbar; colormap('hot');

figure;
imagesc(B1map(:,:,sl)); caxis([0.4 1.5]); colorbar; colormap('hot');




