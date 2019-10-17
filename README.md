# fMRIDec_Python

FastMapping_Tool Version 1.0 (Jun 5, 2018)


Together with this Readme.txt are the folders v_Matlab and v_Python
for Paradigm free Mapping of fMRI data using Lasso and fused Lasso


Before running Matlab codes, please first make sure that SPM (http://www.fil.ion.ucl.ac.uk/spm/) is installed and added to the Matlab path. 

Before running the Python codes, please make sure that Nipy modules (http://nipy.org/nipy/) and the wavelet package PyWavelets (https://github.com/PyWavelets) are downloaded. The FastMapping_Tool module is called FastMap.py

@developper: Wavelets are used only for the noise estimation process. You can use your favorite noise estimation procedure instead.


The folder "DataExample" provides examples of fMRI data used in the main code

The outputs are saved in the OUTPUT folder



%%%%% Run the code on example data by going to the folders:

v_Matlab: matlab main.m       (tested on Matlab 2016b)
v_Python: python main.py      (tested on Python 2.7)

If any problem, please contact Younes Farouj via younes.farouj@epfl.ch.
