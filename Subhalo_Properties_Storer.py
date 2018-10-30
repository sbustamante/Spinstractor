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
#DataFolder = '/home/bustamsn/PhD/DataHaswell/'
DataFolder = '/home/bustamsn/PhD/DataFreya/'
#DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/'
#DataFolder = '/home/bustamsn/PhD/CosmologicalBlackHoles/Data'
DataResults = '.'
Resolution = 512
Simulation = 'cosmobh04/'

nshells = 60
zcut =  1                     # vertical cut in Kpc


#Morphology estimator       # Circularity    Fitting
estimator = 'Circularity'


#========================================================================================
#   GALAXY PLOT
#========================================================================================
attrs = ['pos', 'vel', 'mass', 'age', 'pot', 'id']
sf = load_subfind( 15, dir='%s/Sims%3d/%s/output/'%(DataFolder, Resolution, Simulation), hdf5=True, loadonly=['smty', 'SubhaloCM', 'flty', 'fnsh', 'slty', 'fpos', 'frc2', 'svel', 'spos', 'ffsh','fmty', 'ssph', 'shmt'] )
#Calculating indexes
IDs = np.arange(len(sf.data['fmty']))
mask_bh = sf.data['fmty'][:,5]>0
IDs_sub = np.insert(np.cumsum(sf.data['fnsh'])[:-1], 0, 0)

IDs = IDs[mask_bh]


IDs_sub = IDs_sub[mask_bh]
StellarMass2 = sf.data['fmty'][mask_bh][:,4]
StellarMass = sf.data['smty'][IDs_sub][:,4]
StellarHalfRadius = sf.data['smty'][IDs_sub][:,4]
GasMass = sf.data['smty'][IDs_sub][:,0]
DMMass = sf.data['smty'][IDs_sub][:,1]
Colorgr = sf.data['ssph'][IDs_sub][:,4] - sf.data['ssph'][IDs_sub][:,5]
Spos = sf.data['spos'][IDs_sub]


mask_ord = np.argsort(StellarMass2)[::-1]

IDs = IDs[mask_ord]
IDs_sub = IDs_sub[mask_ord]
StellarMass = StellarMass[mask_ord]
StellarHalfRadius = StellarHalfRadius[mask_ord]
GasMass = GasMass[mask_ord]
DMMass = DMMass[mask_ord]
Colorgr = Colorgr[mask_ord]
StellarMass2 = StellarMass2[mask_ord]
Spos = Spos[mask_ord]


#DT = abs(np.loadtxt( './data/morphologies_cosmobh04.txt' ))
