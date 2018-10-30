#========================================================================================#========================================================================================
#   LIBRARIES
#========================================================================================
import BH_Physics as BHP
from fitool import fit

import numpy as np
from scipy.interpolate import interp1d
from loadmodules import *
from parse_particledata import parse_particledata
import os
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import rc
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.close('all')



#========================================================================================
#   PARAMETERS
#========================================================================================
#Data folder
#DataFolder = '/home/bustamsn/PhD/Data/cosmological_BH/'
#DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/'
DataFolder = '/home/bustamsn/PhD/CosmologicalBlackHoles'
DataResults = '.'
Resolution = 512
Simulation = 'Data'

nshells = 60
zcut =  1                     # vertical cut in Kpc


#========================================================================================
#   GALAXY PLOT
#========================================================================================
attrs = ['pos', 'vel', 'mass', 'age', 'pot', 'id']
s = gadget_readsnap( 15, snappath='%s/%s/snapdir_015'%(DataFolder, Simulation), snapbase='snap_', hdf5=True, loadonlytype=[4,5], loadonly=attrs )
sf = load_subfind( 15, dir='%s/%s/'%(DataFolder, Simulation), hdf5=True, loadonly=['smty', 'flty', 'fnsh', 'slty', 'fpos', 'frc2', 'svel', 'spos', 'ffsh','fmty'] )
#Calculating indexes
IDs = np.arange(len(sf.data['fmty']))
mask_bh = sf.data['fmty'][:,5]>0
IDs = IDs[mask_bh]
StellarMass = sf.data['fmty'][mask_bh][:,4]
mask_ord = np.argsort(StellarMass)[::-1]
StellarMass = StellarMass[mask_ord]
IDs = IDs[mask_ord]


IDHALO = 0
plt.close('all')
s = gadget_readsnap( 15, snappath='%s/%s/snapdir_015'%(DataFolder, Simulation), snapbase='snap_', hdf5=True, loadonlytype=[4,5], loadonly=attrs )
s.calc_sf_indizes( sf )
s.select_halo( sf, use_principal_axis=False, use_cold_gas_spin=False, do_rotation=False, haloid=IDs[IDHALO] )

g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][IDs[IDHALO]] )
g.prep_data()

sdata = g.sgdata['sdata']
mass = sdata['mass']
pos = sdata['pos']
pot = sdata['pot']



s = gadget_readsnap( 15, snappath='%s/%s/snapdir_015'%(DataFolder, Simulation), snapbase='snap_', hdf5=True, loadonlytype=[4,5], loadonly=attrs )
s.calc_sf_indizes( sf )
s.select_halo( sf, use_principal_axis=False, use_cold_gas_spin=False, do_rotation=False, haloid=IDs[IDHALO], subhalo=20 )

g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][IDs[IDHALO]] )
g.prep_data()

sdata = g.sgdata['sdata']
mass1 = sdata['mass']
pos1 = sdata['pos']
pot1 = sdata['pot']




#3D plot

fig = plt.figure()
ax = fig.gca(projection='3d')

mask = g.sgdata['sdata']['eps2'] > 0.7

#ax.plot( pos[mask][:,1], pos[mask][:,2], pos[mask][:,0], 'o', ms=2, color='red' )
#ax.plot( pos[mask==False][:,1], pos[mask==False][:,2], pos[mask==False][:,0], '.', ms=1 )

ax.plot( pos[:,1], pos[:,2], pos[:,0], 'o', ms=2, color='red' )
ax.plot( pos1[:,1], pos1[:,2], pos1[:,0], '.', ms=2, color='blue' )

ax.set_aspect(1)
plt.show()
