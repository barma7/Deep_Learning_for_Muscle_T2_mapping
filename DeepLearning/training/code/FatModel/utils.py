#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:48:51 2018

@author: marco.barbieri21
"""

import numpy as np
from os.path import join
import scipy.io  
import random
# general imports
from config import data_path, B, s, r, exp_folder, Fixed_FF, FF, FFmin, FFmax
from config import Normalization, Nb_pulses, AddNoise, SNRmax, SNRmin
from config import training_size, shuffle

# READ GT-FILE once and for all
GT_ = np.float32(np.loadtxt(join(data_path,exp_folder,'LUT.txt'), delimiter=','))

GT = np.copy(GT_)
GT[:, 0] = GT[:, 0]/100 # T2w, divided by 100 to have values similar to B1 ratio
GT[:, 1] = GT[:, 1]/100 # T2f, divided by 100 to have values similar to B1 ratio

# split data in train and test (with fixed seed)
n = GT.shape[0]
x = np.array(range(n))
if shuffle:
    random.shuffle(x, lambda: s)
    
if training_size == False:
    x_train = x[0:int(r * n)]
    x_test = x[int(r * n):]
else:
    x_train = x[0:training_size]
    x_test = x[training_size:]

## LOAD DATASET from matlab file
Data_w = scipy.io.loadmat(join(data_path,exp_folder, 'dictionary_water.mat'))
Data_w = np.complex64(Data_w['dictionary_water'][:,0:Nb_pulses])
Data_f = scipy.io.loadmat(join(data_path,exp_folder, 'dictionary_fat.mat'))
Data_f = np.complex64(Data_f['dictionary_fat'][:,0:Nb_pulses])

# DATA GENERATOR - queues data for keras.fit_generator
def data_generator(train=True):
    while True:
        input_, labels_ = get_data(train=train)
        yield input_, labels_

# ACTUAL FUNCTION PROVIDING THE DATA
def get_data(train=True):

    # sample B indexes from training/testing subset
    x_ = x_train if train else x_test
    x_sampled = random.sample(list(x_), np.minimum(B, x_.shape[0]))

    # prepare outputs
    data_tensor_w = Data_w[x_sampled]
    data_tensor_f = Data_f[x_sampled]

    labels_tensor = GT[x_sampled]

    ## DATA AGUMENTATION
    if Fixed_FF:
        data_tensor = (1 - FF) * data_tensor_w + FF * data_tensor_f
    else:
        sampled_ff = np.random.uniform(low=FFmin, high=FFmax, size=(len(x_sampled), 1))
        labels_tensor = np.concatenate((sampled_ff, labels_tensor), axis=1)
        data_tensor = (1-sampled_ff)*data_tensor_w + sampled_ff*data_tensor_f

    if AddNoise:
        sampled_snr = np.random.randint(SNRmin, SNRmax+1,len(x_sampled))
        data_tensor = add_gaussian_noise_to_signal(data_tensor, SNR_lin=sampled_snr)

        #sampled_var = np.random.uniform(var_add_min, var_add_max, len(x_sampled))
        #data_tensor = add_gaussian_noise_to_signal_var(data_tensor, var=sampled_var)

    ## DATA PREPROCESSING: STEP WE'RE GOING TO APPLY ALSO TO REAL DATA
    # Normalization
    if Normalization:
        data_tensor = np.divide(data_tensor,np.linalg.norm(data_tensor,axis=1)[:,None])

    data_tensor_ = np.abs(data_tensor)
        
    return data_tensor_, labels_tensor

## Since the singal is a complex signal, we need to make sure half of the variance goes for real and half for imaginary.
def add_gaussian_noise_to_signal(signal, SNR_lin):
    
    L = signal.shape[1]
    n0 = np.divide(np.divide(np.sum(np.square(np.abs(signal)), axis=1), L), SNR_lin)
    
    if np.iscomplexobj(signal):
        n1 = np.multiply(np.random.randn(signal.shape[0], signal.shape[1]), np.sqrt(np.divide(n0, 2))[:, None])
        n2 = np.multiply(np.random.randn(signal.shape[0], signal.shape[1]), np.sqrt(np.divide(n0, 2))[:, None])
        
        noisy_signal = np.complex64(signal + (n1 + 1j*n2))
    else:
        n1 = np.random.randn(signal.shape[0], signal.shape[1]) * np.sqrt(n0)[:, None]
        noisy_signal = np.float32(signal + n1)
        
    return noisy_signal


def add_gaussian_noise_to_signal_var(signal, var):
    L = signal.shape[1]
    n0 = var

    if np.iscomplexobj(signal):
        n1 = np.multiply(np.random.randn(signal.shape[0],signal.shape[1]), np.sqrt(np.divide(n0,2))[:, None])
        n2 = np.multiply(np.random.randn(signal.shape[0],signal.shape[1]), np.sqrt(np.divide(n0,2))[:, None])

        noisy_signal = np.complex64(signal + (n1 + 1j * n2))
    else:
        n1 = np.random.randn(signal.shape[0], signal.shape[1]) * np.sqrt(n0)[:, None]
        noisy_signal = np.float32(signal + n1)

    return noisy_signal

