#==================================================================================================
#       LIBS
#==================================================================================================
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import chain
import numpy as np
from scipy.interpolate import interp1d
import os
import h5py
import arepo
import illustris_python as il
from matplotlib import rc
rc('font',family='Palatino')
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


#==================================================================================================
#       INPUT DATA
#==================================================================================================

#Data folder
#DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/'
DataFolder = '/home/bustamsn/PhD/DataHaswell/'
#DataFolder = '/home/bustamsn/PhD/Data/cosmological_BH/'
#Resolution
Resolution = 512
#Simulation
Simulation = 'cosmobh00'

#==================================================================================================
#       PLOTTING
#==================================================================================================
basePath = '%s/Sims%3d/%s'%(DataFolder, Resolution, Simulation)
#stars = il.snapshot.loadSubhalo(basePath,15,0,'stars')

sub = arepo.Subfind('%s/Sims%3d/%s/output/groups_015'%(DataFolder, Resolution, Simulation), snapshot=15, combineFiles=True )
index = np.argsort(sub.SubhaloMass)[::-1]
Halo = index[15]

filter = [arepo.filter.Halo(sub,subhalo=Halo) ]

sn = arepo.Snapshot('%s/Sims%3d/%s/output/snapdir_015/snap_015.0.hdf5'%(DataFolder, Resolution, Simulation), parttype=[0,5], filter=filter, combineFiles=True, verbose=True)
coor = sn.part0.Coordinates
coorBH = sn.part5.Coordinates

ax = plt.subplot( projection='3d' )
plt.plot( coor[::,0], coor[::,1], coor[::,2], '.', ms=0.1 )
plt.plot( coorBH[::,0], coorBH[::,1], coorBH[::,2], '.', ms=10, color='red' )

#ax.set_xlim( sub.GroupCM[Halo,0]-1000, sub.GroupCM[Halo,0]+1000 )
#ax.set_ylim( sub.GroupCM[Halo,1]-1000, sub.GroupCM[Halo,1]+1000 )
#ax.set_zlim( sub.GroupCM[Halo,2]-1000, sub.GroupCM[Halo,2]+1000 )
ax.set_aspect('equal', 'datalim')    
    
plt.show()
 
 
