
# Import statements
import tensorflow as tf
import math
from keras import backend as K
from keras import optimizers
from keras import metrics
from pathlib import Path
from os.path import join

import h5py

from keras.layers import Input, Dense, LeakyReLU
from keras.models import Model
from keras.callbacks import CSVLogger, ModelCheckpoint
# import general configuration info
from utils import data_generator
from config import N, N_out, lr_init, dev_set, homeDir, save_path


# SETUP AND PRINT NETWORK MODEL
# ####################################################################

def network_model3():

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


input_stream, predictions = network_model3() # create the tensors by the defined function network_model()
model = Model(inputs=(input_stream), outputs=(predictions))  #create the model to be trained by Keras
model.summary()

# DEFINE PROBLEM SPECIFIC LOSS FUNCTION AND COMPILE THE MOEDL
# ####################################################################
Path(join(save_path, "log_loss")).mkdir(parents=True, exist_ok=True)
Path(join(save_path, "weights")).mkdir(parents=True, exist_ok=True)

adam = optimizers.Adam(lr=lr_init, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.0)
model.compile(loss='mae', optimizer=adam, metrics=['mape'])
checkpoint = ModelCheckpoint(join(save_path, 'weights', 'MODEL.{epoch:02d}.hdf5'), save_weights_only=True, period=250)

# SETUP AND COMPILE NETWORK MODEL
# ####################################################################
#lr_scheduler_cb = LearningRateScheduler(lambda epoch: lr_init)
csv_logger = CSVLogger(join(save_path, 'log_loss', 'log_loss_MODEL.csv'), separator=',', append=False)

if dev_set:
    model.fit_generator(data_generator(train=True), steps_per_epoch=1000, epochs=500, verbose=1, callbacks=[csv_logger, checkpoint],
                        validation_data=data_generator(train=False), validation_steps=500, initial_epoch=0)
else:
    model.fit_generator(data_generator(train=True), steps_per_epoch=1000, epochs=500, verbose=1, callbacks=[csv_logger, checkpoint],initial_epoch=0)

