#========================================================================================
#   LIBRARIES
#========================================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
import time
plt.close('all')


#========================================================================================
#   PARAMETERS
#========================================================================================
#Data folder
DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/Sims256/'
#Simulation
Simulation = 'cosmobh05'
#Number of chunks (same number of used processors)
N_proc = 256

#========================================================================================
#   Extracting data
#========================================================================================
indexes = np.loadtxt('%s%s/analysis/BH_IDs.txt'%(DataFolder,Simulation))[:,[0,1]]

os.system('rm tmp.txt')
for i in xrange(N_proc):
    print 'In file', i
    #Loading data
    str_cmd = "less %s%s/output/blackhole_mergers/blackhole_mergers_%d.txt >> tmp.txt"%(DataFolder,Simulation,i)
    os.system(str_cmd)

#Storing concatenated file    
data = np.loadtxt('tmp.txt')
os.system('rm tmp.txt')
f= open('%s%s/analysis/BH_Mergers.txt'%(DataFolder,Simulation), 'a')
np.savetxt(f, data[:,1:], fmt='%e %d %e %d %e')
f.close() 
 
