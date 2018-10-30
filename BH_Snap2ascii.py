#==================================================================================================
#       LIBS
#==================================================================================================
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import chain
import numpy as np
from scipy.interpolate import interp1d
import os
import h5py
import arepo 
#==================================================================================================
#       INPUT DATA
#==================================================================================================

#Data folder
#DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/'
DataFolder = '/home/bustamsn/PhD/DataHaswell/'
#DataFolder = '/home/bustamsn/PhD/DataFreya/'
#Resolution
Resolution = 512
#Simulation
#Simulations = ['cosmobh00', 'cosmobh01', 'cosmobh02', 'cosmobh03', 'cosmobh04']
Simulations = ['cosmobh02']

#Total bins
total_bins = 20
#Interpolation type
interptype = 'linear'

#Observational Kormendi-Ho function
Mbulg = np.linspace( 8, 12, 100 )
Mbh = 1.06352*( Mbulg - np.log10(2.7e8) ) + 6

#Colors
colors = ['black', 'blue', 'red', 'green', 'orange', 'gray', 'cyan']

#==================================================================================================
#       PLOTTING
#==================================================================================================
plt.close('all')
plt.figure( figsize = (9,14) )
    
plt.subplot( 5, 2, 1 )    

for i_sim, sim in enumerate(Simulations):

    for snapshot in [8,]:
        
        #Loading subfind halos
        sub = arepo.Snapshot( '%s/Sims%3d/%s/output/snapdir_015'%(DataFolder, Resolution, sim), snapshot=snapshot, combineFiles=True, parttype=[5] )
        
        try:
            os.system('mkdir %s/Sims%3d/%s/analysis/BH_snap_%d'%(DataFolder,Resolution,sim,snapshot))
        except:
            print 'Not possible to create folder %s/Sims%3d/%s/analysis/BH_snap_%d'%(DataFolder,Resolution,sim,snapshot)
        
        for var in dir(sub.part5):
            try:
                np.savetxt( '%s/Sims%3d/%s/analysis/BH_snap_%d/%s.txt'%(DataFolder,Resolution,sim,snapshot, var), sub.part5[var] )
            except:
                None
        
        
        
