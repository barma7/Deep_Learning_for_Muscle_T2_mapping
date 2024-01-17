# Deep_Learning_for_Muscle_T2_mapping
The repository provides the source code for "A Deep Learning Approach for Fast Muscle Water T2 Mapping with Subject Specific Fat T2 Calibration from Multi-Spin-Echo Acquisitions"

## Deep Learning Application ##
The deep learning application has been developed in TensorFlow.
TensorFlow, Numpy, and Scipy are required to run the deep learning application.

The EPG code for simulating the Multi-Echo-Spin-Echo (MESE) data is provided in Matlab.
The core EPG function is taken from the StimFit toolbox (https://github.com/usc-mrel/StimFit).
Since the actual pulse shapes used by the vendors in implementing the MESE sequence are proprietary information, the simulations are provided for SINC pulses with TBW = 2.

The functions used to generate the excitation and refocusing pulses and respective slice profiles are provided and make use of the RF Pulse Design toolbox created by John Pauly (https://rsl.stanford.edu/research/software.html) 

## DATA ##
The raw MESE data (in NIfTI format) and processed data can be downloaded from Zenodo (https://doi.org/10.5281/zenodo.10520542) 

## Alternative T2 fitting methods ##
Conventional EPG fitting algorithms based on Dictionary and Non-Linear Lease-Squared (NLSQ) methods are also provided in Matlab.




