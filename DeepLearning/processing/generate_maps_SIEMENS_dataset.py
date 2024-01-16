
import numpy as np
from keras.layers import Input, Dense, Concatenate, LeakyReLU, Dropout
from keras.models import Model
import dosma as dm
import os
import matplotlib.pyplot as plt
from pathlib import Path
import time
import pandas as pd
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from scipy.stats import gmean
from os.path import join

# LOAD AND PREPARE DATA

# NN training folders
homeDir = r'C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PaperSharing'
dataset_folder = join(homeDir, "DATA", "SIEMENS")
subjects = ["02", "03", "04", "05"]

ffmap_names = ["ff_s1.nii.gz"]
t2map_names = ["t2_s1.nii.gz"]
b1map_names = ["b1_s1.nii.gz"]
tse_names = ["mese_vol.nii.gz"]

sv_root = "MAPS/DL/with_b1/FatNet4wb1-snr_1e3-5e4"
sv_fldr = "MuscleNet"
sv_name = join(sv_root, sv_fldr)

## NN Parameters
# model3 W/ B1 estimation
wieghts_path = join(homeDir, r"DeepLearning\training\results\SIEMENS\ETL_17_ESP_7.5\MuscleNet\Noise\var5e-5_1e-9\t2w-10-110-t2f-50-250-b1-0.3-1.2")
# FatNet W/ B1 estimation
FatModel_weights =join(homeDir, r"DeepLearning\training\results\SIEMENS\ETL_17_ESP_7.5\FatNet/Noise/snr-1e3_5e4\t2w-10-110-t2f-50-250-b1-0.3-1.2")

# Preprocessing Steps
Normalization = True
b1_estimation = True
b1_estimation_fat = True

# Sequence parameters
ETL = 17 # Echo Train Length

nr = dm.NiftiReader()
nwr = dm.NiftiWriter()

## DEFINE NN MODEL For T2water and FF and LOAD WEIGHTS
Nb_pulses = 17

# MuslceNet params
if b1_estimation:
    N_out = 3
else:
    N_out = 2

N = Nb_pulses + 1  # input size

def muscle_net3():

    my_input = Input(shape=(N,),name='input')
    x = Dense(512, kernel_initializer='glorot_normal', name='hidden0')(my_input)
    x = LeakyReLU(name='hidden0_1')(x)
    x = Dense(256, kernel_initializer='glorot_normal', name='hidden1')(x)
    x = LeakyReLU(name='hidden1_1')(x)
    x = Dense(128, kernel_initializer='glorot_normal', name='hidden2')(x)
    x = LeakyReLU(name='hidden2_1')(x)
    x = Dense(64, kernel_initializer='glorot_normal', name='hidden3')(x)
    x = LeakyReLU(name='hidden3_1')(x)
    x = Dense(N_out, activation='linear', kernel_initializer='glorot_normal', name='output')(x)

    my_output = x

    return my_input, my_output

input_stream, predictions = muscle_net3()  # create the tensors by the defined function network_model()
MuscleNet = Model(inputs=(input_stream), outputs=(predictions))  # create the model to be trained by Keras
MuscleNet.load_weights(join(wieghts_path, 'weights', 'MODEL.250.hdf5'))

# FatNet parameters

# MuslceNet params
if b1_estimation_fat:
    N_out_fat = 3
else:
    N_out_fat = 2

N_fat = Nb_pulses # input size

def FatNet3():

    my_input = Input(shape=(N_fat,),name='input')
    x = Dense(512, kernel_initializer='glorot_normal', name='hidden0')(my_input)
    x = LeakyReLU(name='hidden0_1')(x)
    x = Dense(256, kernel_initializer='glorot_normal', name='hidden1')(x)
    x = LeakyReLU(name='hidden1_1')(x)
    x = Dense(128, kernel_initializer='glorot_normal', name='hidden2')(x)
    x = LeakyReLU(name='hidden2_1')(x)
    x = Dense(64, kernel_initializer='glorot_normal', name='hidden3')(x)
    x = LeakyReLU(name='hidden3_1')(x)
    x = Dense(N_out_fat, activation='linear', kernel_initializer='glorot_normal', name='output')(x)

    my_output = x

    return my_input, my_output

input_stream_, predictions_ = FatNet3()  # create the tensors by the defined function network_model()
FatNet = Model(inputs=(input_stream_), outputs=(predictions_))
FatNet.load_weights(join(FatModel_weights, 'weights', 'MODEL.500.hdf5'))

strt_tot = time.time()
##  Reconstruct volume from dicom and load medical volume for saving nifti files
for sub in subjects:
    print('processing subject ',sub)
    for stack in range(0, len(tse_names)):
        print('processing stack ', stack)

        nn_svfldr = join(dataset_folder, sub, sv_root)
        SvFldr = join(dataset_folder, sub, sv_name)
        Path(SvFldr).mkdir(parents=True, exist_ok=True)

        #READ DATA
        stack_nifti = nr.load(join(dataset_folder, sub, tse_names[stack]))
        Data = stack_nifti.volume
        
        # extract dimension of the data
        nb_row = Data.shape[0]
        nb_col = Data.shape[1]
        nb_slices = Data.shape[2]
        nb_echoes = Data.shape[3]

        #check if T2fat exist in folder
        if not os.path.exists(join(nn_svfldr, "Fat_calibration.csv")):
            print('No fat calibration in folder ', nn_svfldr, ' Processing Fat Calibration...')
            fat_nifti = nr.load(join(dataset_folder, sub, tse_names[0]))
            fat_mask = nr.load(join(dataset_folder, sub, "fat_mask.nii"))
            
            # extract dimension of the data
            nb_row = fat_nifti.volume.shape[0]
            nb_col = fat_nifti.volume.shape[1]
            nb_slices = fat_nifti.volume.shape[2]
            nb_echoes = fat_nifti.volume.shape[3]

            fat_Data = fat_nifti.volume
            
            # flat sub_fat_data
            fat_data_flat = np.reshape(fat_Data, (nb_row*nb_col*fat_Data.shape[2], nb_echoes))
            fat_mask_flat = fat_mask.volume.flatten()

            fat_data = fat_data_flat[fat_mask_flat > 0]

            if Normalization:
                fat_data = np.divide(fat_data, np.linalg.norm(fat_data, axis=1)[:, None])

            strt = time.time()
            cal_prediction = FatNet.predict(fat_data)
            stp = time.time()
            print(stp - strt)

            pred = np.array(cal_prediction)

            #pred[:, 0] = pred[:, 0]*1000
            pred[:, 0] = pred[:, 0]*100
            pred[:, 1] = pred[:, 1]*100
            
            # cropping values
            pred[pred < 0] = 0
            pred[pred[:, 1] > 250, 1] = 250
            pred[pred[:, 0] > 110, 0] = 110
            pred[pred[:, 1] < 50, 1] = 50
            pred[pred[:, 0] < 10, 0] = 10
            pred[pred[:, 2] > 1.2, 2] = 1.2
            pred[pred[:, 2] < 0.3, 2] = 0.3
            
            pred = pred[~np.isnan(pred).any(axis=1), :]

            t2_water_scatter = np.copy(pred[:,0])
            t2_fat_scatter = np.copy(pred[:, 1])
            b1_fat_scatter = np.copy(pred[:, 2])

            # get t2_water and t2_fat values where b1_fat > 0.7
            t2_water_scatter = t2_water_scatter[b1_fat_scatter > 0.7]
            t2_fat_scatter = t2_fat_scatter[b1_fat_scatter > 0.7]
            
            t2_water_scatter = t2_water_scatter.reshape(-1, 1)
            t2_fat_scatter = t2_fat_scatter.reshape(-1, 1)
            
            # Create scatter data structure
            t2_predicted = np.concatenate((t2_water_scatter, t2_fat_scatter), axis=1)
            t2_fat_mean = np.nanmean(t2_fat_scatter)
            t2_water_mean = np.nanmean(t2_water_scatter)

            t2_fat_geom = gmean(t2_fat_scatter, nan_policy='omit')
            t2_water_geom = gmean(t2_water_scatter, nan_policy='omit')

            # np.save(os.path.join(nn_svfldr, "scatter_t2w_t2f.npy"), t2_predicted)

            T2fat = t2_fat_mean
            
            fat_dict = {'T2water_mean': [t2_water_mean], 'T2fat_mean': [t2_fat_mean],
                        'T2water_geom': t2_water_geom, 'T2fat_geom': t2_fat_geom}

            fat_dt = pd.DataFrame(fat_dict)
            fat_dt.to_csv(os.path.join(nn_svfldr, "Fat_calibration.csv"))

            print("T2 fat mean: ", t2_fat_mean, " T2 water mean: ", t2_water_mean,
                  " T2 fat geom:", t2_fat_geom, " T2 water geom: ", t2_water_geom)

            # ## FIND PDF OF T2fat and T2water
            # print('Creatuing pdfs...')
            # grid_fat = GridSearchCV(KernelDensity(),
            #                         {'bandwidth': np.linspace(0.1, 5.0, 70)},
            #                         cv=20, n_jobs=8)  # 20-fold cross-validation
            # grid_water = GridSearchCV(KernelDensity(),
            #                           {'bandwidth': np.linspace(0.1, 5.0, 70)},
            #                           cv=20, n_jobs=8)  # 20-fold cross-validation

            # grid_fat.fit(pred[:, 1].reshape(-1, 1))
            # grid_water.fit(pred[:, 0].reshape(-1, 1))

            # print(grid_fat.best_params_, grid_water.best_params_)

            # kde_fat = KernelDensity(kernel='gaussian', bandwidth=grid_fat.best_params_.get('bandwidth')).fit(
            #     pred[:, 1].reshape(-1, 1))
            # kde_water = KernelDensity(kernel='gaussian', bandwidth=grid_water.best_params_.get('bandwidth')).fit(
            #     pred[:, 0].reshape(-1, 1))

            # xx_fat = np.linspace(50,200,150).reshape(-1,1)
            # xx_water = np.linspace(0,80,150).reshape(-1,1)

            # logprob_f = kde_fat.score_samples(xx_fat)
            # logprob_w = kde_water.score_samples(xx_water)

            # pdf_f = np.exp(logprob_f)
            # pdf_w = np.exp(logprob_w)

            # dist_fat = np.concatenate((xx_fat, pdf_f.reshape(-1,1)), axis=1)
            # dist_water = np.concatenate((xx_water, pdf_w.reshape(-1,1)), axis=1)

            # np.save(os.path.join(nn_svfldr, "dist_t2fat.npy"), dist_fat)
            # np.save(os.path.join(nn_svfldr, "dist_t2water_fat.npy"), dist_water)

            # T2fat_kde = xx_fat[np.argmax(pdf_f)][0]
            # T2water_kde = xx_water[np.argmax(pdf_w)][0]

            # fsize = 34
            # params = {'legend.fontsize': 28,
            #           'axes.labelsize': fsize,
            #           'axes.titlesize':fsize,
            #           'xtick.labelsize':fsize,
            #           'ytick.labelsize':fsize}
            # plt.rcParams.update(params)
            # plt.rc('text', usetex=True)
            # plt.rc('axes', linewidth=2)
            # plt.rc('font', weight='bold')

            # plt.figure(figsize=(15,8))
            # plt.plot(xx_fat, pdf_f, linewidth=4)
            # plt.xlabel(r'\textbf{$T_2$ (ms)}')
            # plt.ylabel(r'\textbf{Probability Density}')
            # plt.savefig(os.path.join(nn_svfldr,'T2fat_distribution.png'), format='png', dpi=300)
            # #plt.show()

            # plt.figure(figsize=(15,8))
            # plt.plot(xx_water, pdf_w, linewidth=4)
            # plt.xlabel(r'\textbf{$T_2$ (ms)}')
            # plt.ylabel(r'\textbf{Probability Density}')
            # plt.savefig(os.path.join(nn_svfldr,'T2water_fat_distribution.png'), format='png', dpi=300)
            # #plt.show()
        else:
            print('Found fat calibration in folder ', nn_svfldr, ' Processing Fat Calibration skipped')
            fat_calibration = pd.read_csv(os.path.join(nn_svfldr, "Fat_calibration.csv"))
            T2fat = float(fat_calibration["T2fat_mean"][0])


        ## QUANTIATIVE MAPS PROCESSING
        print('Computing quantiative maps')

        FFmap = np.reshape(np.zeros((nb_row, nb_col, nb_slices)), (nb_row*nb_col*nb_slices))
        T2map = np.reshape(np.zeros((nb_row, nb_col, nb_slices)), (nb_row*nb_col*nb_slices))
        if b1_estimation:
            B1map = np.reshape(np.zeros((nb_row, nb_col, nb_slices)), (nb_row*nb_col*nb_slices))

        strt1 = time.time()

        data_temp = np.copy(np.reshape(Data, (nb_row*nb_col*nb_slices, ETL)))
        non_zero_indexes = np.where(data_temp[:, 0] > 50)[0]
        data_temp = data_temp[non_zero_indexes, :]

        t2fat = T2fat/100 * np.ones((data_temp.shape[0], 1))

        if Normalization:
            data_temp = np.divide(data_temp, np.linalg.norm(data_temp, axis=1)[:, None])

        data_temp = np.concatenate((t2fat, data_temp), axis=1)

        strt = time.time()
        prediction = MuscleNet.predict(data_temp)
        stp = time.time()
        #print(stp - strt)

        pred = np.array(prediction)
        # cropping values less than zero
        pred[pred < 0] = 0
        pred[:, 1] = pred[:, 1]*100

        # cropping values to correct ranges
        pred[np.where(pred[:, 1] > 80), 1] = 80
        pred[np.where(pred[:, 1] < 10), 1] = 10
        pred[np.where(pred[:, 0] > 1), 0] = 1
        pred[np.where(pred[:, 0] < 0), 0] = 0
        if b1_estimation:
            pred[np.where(pred[:, 2] > 1.2), 1] = 1.2
            pred[np.where(pred[:, 2] < 0.3), 1] = 0.3

        FFmap[non_zero_indexes] = pred[:, 0]
        T2map[non_zero_indexes] = pred[:, 1]
        if b1_estimation:
            B1map[non_zero_indexes] = pred[:, 2]

        FFmap = np.reshape(FFmap, (nb_row, nb_col, nb_slices))  
        T2map = np.reshape(T2map, (nb_row, nb_col, nb_slices))
        if b1_estimation:
            B1map = np.reshape(B1map, (nb_row, nb_col, nb_slices))

        print("All slices processed in ", stp - strt)
        print("writing volumes to nifti ")
        FFmap_mv = dm.MedicalVolume(100*FFmap, stack_nifti.affine)
        T2map_mv = dm.MedicalVolume(T2map, stack_nifti.affine)
        if b1_estimation:
            B1map_mv = dm.MedicalVolume(B1map, stack_nifti.affine)

        nwr.save(FFmap_mv, os.path.join(SvFldr, ffmap_names[stack]))
        nwr.save(T2map_mv, os.path.join(SvFldr, t2map_names[stack]))
        if b1_estimation:
            nwr.save(B1map_mv, os.path.join(SvFldr, b1map_names[stack]))

        del stack_nifti

stop_tot = time.time()
print("Dataset processed in: ", stop_tot - strt_tot)





