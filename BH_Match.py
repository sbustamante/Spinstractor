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
import arepo
plt.close('all')


#========================================================================================
#   PARAMETERS
#========================================================================================
#Resolution
Resolution = 512
#Data folder
#DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/'
#DataFolder = '/home/bustamsn/bustamsn_haswell/Sims%d/'%(Resolution)
DataFolder = '/home/bustamsn/PhD/DataHaswell/Sims%d/'%(Resolution)
#Simulation Reference
Sim_ref = 'cosmobh00'
#Simulation
Simulations = ['cosmobh01', 'cosmobh02']
#Number of chunks (same number of used processors)
N_proc = 512

#========================================================================================
#   Extracting data
#========================================================================================

indexes = np.loadtxt('%s%s/analysis/BH_IDs.txt'%(DataFolder,Sim_ref))[:,[0,1,2,3]]
#Detecting BHs at z=0
mask = indexes[:,2] == 1.0
indexes = indexes[mask]
mask_sort = np.argsort(indexes[:,1])
indexes = indexes[mask_sort]

coords = np.loadtxt('%s%s/analysis/BH_coord.txt'%(DataFolder,Sim_ref))
mask_evolution = np.in1d(coords[:,0].astype(int),indexes[:,1].astype(int))
coords = coords[mask_evolution]
mask_sort = np.argsort(coords[:,0])
coords = coords[mask_sort]

ID_BHS = []
ID_SIM = []
for Simulation in Simulations:
    ID_bhs = []
    ID_sim = []
    
    #Loading data from simulation to compare
    indexes_c = np.loadtxt('%s%s/analysis/BH_IDs.txt'%(DataFolder,Simulation))[:,[0,1,2]]
    #Detecting BHs at z=0
    mask = indexes_c[:,2] == 1.0
    indexes_c = indexes_c[mask]
    
    coords_c = np.loadtxt('%s%s/analysis/BH_coord.txt'%(DataFolder,Simulation))
    mask_evolution = np.in1d(coords_c[:,0].astype(int),indexes_c[:,1].astype(int))
    coords_c = coords_c[mask_evolution]
    mask_sort = np.argsort(coords_c[:,0])
    coords_c = coords_c[mask_sort]
        
    
    for i in xrange(len(indexes)):
        id_closest = np.argsort(np.linalg.norm( coords[i,[1,2,3]] - coords_c[:,[1,2,3]], axis = 1 ))[0]
        try:
            id_indexes = np.arange(len(indexes_c))[ indexes_c[:,1]==coords_c[id_closest,0] ][0]
            ID_bhs.append( indexes_c[id_indexes,0] )
            ID_sim.append( indexes_c[id_indexes,1] )
        except:
            ID_bhs.append( -1 )
            ID_sim.append( -1 )
        
    ID_BHS.append( ID_bhs )
    ID_SIM.append( ID_sim )
    
ID_BHS = np.array(ID_BHS)
ID_SIM = np.array(ID_SIM)
#Masking with mass
mask_mass = np.argsort(indexes[:,3])[::-1]
np.savetxt( 'Table_comparison.txt', np.array( [ indexes[:,0][mask_mass], ID_BHS[0][mask_mass], ID_BHS[1][mask_mass], \
                                                indexes[:,1][mask_mass], ID_SIM[0][mask_mass], ID_SIM[1][mask_mass] ] ).T, fmt = '%d' )
        
        
 
