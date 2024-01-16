#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:49:41 2018

@author: marco.barbieri21
"""
# DATA PATHS AND PARAMETERS
from os.path import join

homeDir = r'C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PaperSharing\DeepLearning'
save_path = join(homeDir, r'training\results\PHILIPS\ETL_17_ESP_7.6\FatNet\Noise\fixed-FF\snr-e3-e4-t2w-10-110-t2f-50-200')
data_path = r'C:\Users\mb7\OneDrive_Stanford\Research\WorkHome\Projects\Muscle-T2-Mapping\PaperSharing\sim-data\PHILIPS'
exp_folder = r'cpmg_ET_7.6ms_ETL_17\propietary_pulses\100000-examples_t2w-10-110-t2f-50-200-b1-1.2'

# decide to shuffle examples or not
shuffle = False

# if True make the norm of the vector equals to one
Normalization = True

## DEFINE DATA AUGMENTATION STEP
AddNoise = True
SNRmin = 1*10**3
SNRmax = 1*10**4

Fixed_FF = True
FF = 0.9
FFmin = 0
FFmax = 1

## DEFINE INPUT AND OUTPUT DIMENSIONS
# number of MRF pulses
Nb_pulses = 17

Input_channels = 1 # 1 or 2

## NETWORK PARAMS
# copmpute output dimension
N = Nb_pulses
if Fixed_FF:
    N_out = 3
else:
    N_out = 4

## TRAINING PARAMS
# if True test examples are used to estimate validation loss during training
dev_set = True

# if False training/test ratio is defined by the parameter r, otherwise define a number
training_size = False

s = 0.12345  # random seed to ensure reproducibility
r = 9.0/10.0  # training on testing data ratio
B = 500  # batch size
lr_init = 0.0001 # learning rate
