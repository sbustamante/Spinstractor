#========================================================================================
#   LIBRARIES
#========================================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
#import arepo 
import time
plt.close('all')


#========================================================================================
#   PARAMETERS
#========================================================================================
#Data folder
DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/'
#Simulation
Simulation = 'cosmobh05'
#Number of chunks (same number of used processors)
N_proc = 256

#========================================================================================
#   Extracting data
#========================================================================================
indexes = np.loadtxt('%s%s/analysis/BH_IDs.txt'%(DataFolder,Simulation))[:,[0,1]]

os.system('mv %s%s/analysis/spins/* %s%s/analysis/spinsbk/'%(DataFolder,Simulation,DataFolder,Simulation))
for i in xrange(N_proc):
    print 'In file', i
    #Loading data
    str_cmd = "sed 's/^[^=]*=/=/' %s%s/output/blackhole_details/blackhole_spin_%d.txt | cut -c 2- > tmp.txt"%(DataFolder,Simulation,i)
    os.system(str_cmd)
    data = np.loadtxt('tmp.txt')
    os.system('rm tmp.txt')
    
    for i_fl, i_bh in zip(indexes[:,0].astype(int),indexes[:,1].astype(int)):
        #Saving current BH
        data_bh = data[data[:,0] == i_bh]
        data_bh = data_bh[ np.unique(data_bh[:,1],return_index=True)[1] ]
        if data_bh.shape[0] > 0:
            f= open('%s%s/analysis/spins/BHt_%d.txt'%(DataFolder,Simulation,i_fl), 'a')
            np.savetxt(f, data_bh[:,1:], fmt='%e %e %e %e %e %e %e %e %e %e %d')
            f.close()

#Sortering data
for i_bh in indexes[:,0].astype(int):
    try:
        data = np.loadtxt('%s%s/analysis/spins/BHt_%d.txt'%(DataFolder,Simulation,i_bh))
        data = data[np.argsort(data[:,0])]
        np.savetxt('%s%s/analysis/spins/BH_%d.txt'%(DataFolder,Simulation,i_bh), data, fmt='%e %e %e %e %e %e %e %e %e %e %d')
    except:
        np.savetxt('%s%s/analysis/spins/BH_%d.txt'%(DataFolder,Simulation,i_bh), data)
 
