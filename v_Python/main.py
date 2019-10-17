import sys
import numpy as np
import os
# import Fast Mapping module
import FastMap


dir_name = 'DataExample'
base_filename = 'TCN_Visual_Voxel' 
format = 'npy'

TCN = np.load(os.path.join(dir_name, base_filename + "." + format))


dt = 2.04       #  TR
du = 32         #  HRF duration   
T  = TCN.shape[0]

H, Nh = FastMap.MakeHrfToeplitz(T,dt,du)


Nit = 100      # FISTA iterations
condition =  "blocks"       # Blocks or Spikes prior
lam =  1                   # Temporal regularization weight




# Main FISTA loop
TC = FastMap.volumap(TCN, H, Nh, lam, Nit, condition)
    


output_dir = os.path.join('./OUTPUT',base_filename)


if not os.path.exists(output_dir): 
	os.makedirs(output_dir)

np.save(os.path.join(output_dir,'TC.npy') ,TC)	


import pylab as pl
import matplotlib.pyplot as plt

if (len(TCN.shape) == 1):
	plt.plot(TC)
	plt.xlabel('timepoints (seconds / TR)')
	plt.ylabel('activation')
elif (len(TCN.shape) == 2):
    plt.imshow(TC.T,aspect='auto')
    plt.xlabel('timepoints (seconds / TR)')
    plt.ylabel('voxels')

# plt.show()    
pl.savefig(os.path.join(output_dir,'TC.pdf'))



