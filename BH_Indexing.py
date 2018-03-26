#========================================================================================
#   LIBRARIES
#========================================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
plt.close('all')

#========================================================================================
#   PARAMETERS
#========================================================================================
#Data folder
DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/Sims512/'
#Simulation
Simulation = 'cosmobh05'
#Number of chunks (same number of used processors)
N_proc = 256
#Concatenate files
conc = True

#========================================================================================
#   Extracting data
#========================================================================================
#Creating (if not existent) folder with postprocessed data
os.system('mkdir %s%s/analysis'%(DataFolder,Simulation))

if conc:
    os.system('rm %s%s/analysis/massdistro.txt'%(DataFolder,Simulation))
    for i in xrange(N_proc):
        print 'In file', i
        #Loading data
        str_cmd = "sed 's/^[^=]*=/=/' %s%s/output/blackhole_details/blackhole_spin_%d.txt | cut -c 2- > tmp.txt"%(DataFolder,Simulation,i)
        os.system(str_cmd)
        data = np.loadtxt('tmp.txt')
        os.system('rm tmp.txt')
        
        #Saving last times
        data = data[np.argsort(data[:,1])[::-1]]
        data = data[ np.unique(data[:,0],return_index=True)[1] ]
        f= open('%s%s/analysis/massdistro.txt'%(DataFolder,Simulation), 'a')
        np.savetxt(f, data[:,[0,1,2,3,-1]], fmt='%d %e %e %e %d')
    f.close()

data = np.loadtxt('%s%s/analysis/massdistro.txt'%(DataFolder,Simulation))
data = data[np.argsort(data[:,1])[::-1]]
data = data[ np.unique(data[:,0],return_index=True)[1] ]
y = np.zeros([data.shape[0], 6])
y[:,0] = np.arange(data.shape[0])
y[:,1:] = data[np.argsort(data[:,3])[::-1]][:,[0,1,3,2,-1]]
np.savetxt('%s%s/analysis/BH_IDs.txt'%(DataFolder,Simulation), y, fmt='%d %d %f %e %f %d')
 
