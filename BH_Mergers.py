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
Simulation = 'cosmobh03'
#Number of chunks (same number of used processors)
N_proc = 256

#========================================================================================
#   Extracting data
#========================================================================================
indexes = np.loadtxt('%s%s/analysis/BH_IDs.txt'%(DataFolder,Simulation))[:,[0,1]]

#os.system('rm tmp.txt')
for i in xrange(N_proc):
    print 'In file', i
    #Loading data
    str_cmd = "less %s%s/output/blackhole_mergers/blackhole_mergers_%d.txt >> tmp.txt"%(DataFolder,Simulation,i)
    os.system(str_cmd)

#Storing concatenated file    
data = np.loadtxt('tmp.txt')
os.system('rm tmp.txt')
f= open('%s%s/analysis/BH_Mergers_R.txt'%(DataFolder,Simulation), 'a')
np.savetxt(f, data[:,1:], fmt='%e %d %e %d %e')
f.close() 
 
#========================================================================================
#   Applying correction to merger file
#========================================================================================
data_merger = np.loadtxt('%s%s/analysis/BH_Mergers_R.txt'%(DataFolder,Simulation))

M1 = []
for i in xrange(len(data_merger)):
    time = data_merger[i,0]
    id_mr = indexes[indexes[:,1]==data_merger[i,1].astype(int),0].astype(int)
    if len(id_mr)==1:
        try:
            id_mr = id_mr[0]
            #Loading data of current BH
            data_BHi = np.loadtxt('%s%s/analysis/spins/BH_%d.txt'%(DataFolder,Simulation,id_mr))
            mask_t = data_BHi[:,0] > time
            M1.append( data_BHi[mask_t,2][0] - data_merger[i,4] )
        except:
            print i
            M1.append( data_merger[i,2] )
    else:
        M1.append( data_merger[i,2] )
M1 = np.array(M1)

np.savetxt( '%s%s/analysis/BH_Mergers.txt'%(DataFolder,Simulation), np.array( [data_merger[:,0], data_merger[:,1], M1, data_merger[:,3], data_merger[:,4]] ).T, fmt='%e %d %e %d %e' )
        

