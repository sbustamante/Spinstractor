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
from collections import Counter

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import rc
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.close('all')


#========================================================================================
#   FUNCTIONS
#========================================================================================
def sersic_prof1( x, ba, b, n):
    return ba * np.exp(-(x / b)**n )  # warning!! b is not Reff here!!!

def exp_prof(x, da, h):
    return da * np.exp(-x / h ) 

def total_profile(x, da, h, ba, b, n):
    y = exp_prof(x, da, h) + sersic_prof1(x, ba, b, n) 
    return (y)



def Exponential_Disk( r, sigma0, Rd ):
    return sigma0*np.exp( -r/Rd )

def Sersic_Profile( r, sigma_eff, Reff, n ):
    
    b = lambda n: 2*n - 1/3. + 4/(405.*n) + 46/(25515.*n**2) + 131/(1148175*n**3) - 2194697/(30690717750*n**4)
    return sigma_eff*np.exp( -b(n)*( (r/Reff)**(1./n) - 1 ) )

def Total_Profile( r, a, b, c, d, e ):
    return Exponential_Disk( r, a, b ) + Sersic_Profile( r, c, d, e )



#========================================================================================
#   PARAMETERS
#========================================================================================
Simulation0 = 'cosmobh00'
Simulation1 = 'cosmobh01'

#========================================================================================
#   PLOTS
#========================================================================================
StellarMass0 = np.loadtxt( './data/stellarmass_%s.txt'%Simulation0 )
StellarMass1 = np.loadtxt( './data/stellarmass_%s.txt'%Simulation1 )

#GasMass0 = np.loadtxt( './data/gasmass_%s.txt'%Simulation0 )
#GasMass1 = np.loadtxt( './data/gasmass_%s.txt'%Simulation1 )

#DMMass0 = np.loadtxt( './data/dmmass_%s.txt'%Simulation0 )
#DMMass1 = np.loadtxt( './data/dmmass_%s.txt'%Simulation1 )

DT0 = abs(np.loadtxt( './data/morphologies_%s.txt'%Simulation0 ))
DT1 = abs(np.loadtxt( './data/morphologies_%s.txt'%Simulation1 ))

Color0 = np.loadtxt( './data/color_%s.txt'%Simulation0 )
Color1 = np.loadtxt( './data/color_%s.txt'%Simulation1 )

Spos0 = np.loadtxt( './data/spos_%s.txt'%Simulation0 )
Spos1 = np.loadtxt( './data/spos_%s.txt'%Simulation1 )

#Matching
index = np.zeros( len(StellarMass0) )
residual = np.zeros( len(StellarMass0) )
ID1 = np.arange( len(StellarMass1) )
for i,r in enumerate(Spos0):
    index[i] = int(ID1[np.argsort( np.linalg.norm( Spos1 - r, axis=1 ) )[0]])
    residual[i] = np.sort( np.linalg.norm( Spos1 - r, axis=1 ) )[0]
    
index = index.astype(int)
    
#Discarding not-matching haloes
mask_nm = (residual<=1e2)#*(np.log10(StellarMass1[index])+10>9.5)#*(Color0>0.6)


plt.figure( figsize=(15,10) )
mask_nm = (residual<=1e2)*(np.log10(StellarMass1[index])+10>9.5)*(Color1[index]>0.6)
plt.scatter( np.log10(StellarMass0[mask_nm]) + 10 , Color0[mask_nm], alpha = 0.5, marker = 's', s=40, vmin=0, vmax=1, zorder = 10, color='blue')#, c=DT0[mask_nm], cmap = 'jet_r' )
plt.scatter( np.log10(StellarMass1[index[mask_nm]]) + 10 , Color1[index[mask_nm]], alpha = 0.5, marker = 'o', s=40, vmin=0, vmax=1, zorder = 10, color='red')#, c=DT1[index[mask_nm]], cmap = 'jet_r' )
#cb1 = plt.colorbar()
#cb1.set_label( 'D / T', fontsize = 18 )

for i in xrange( sum(mask_nm) ):
    #plt.arrow( [np.log10(StellarMass0[mask_nm])[i] + 10, np.log10(StellarMass1[index[mask_nm]])[i] + 10], [Color0[mask_nm][i], Color1[index[mask_nm]][i]]
    plt.arrow( np.log10(StellarMass0[mask_nm])[i] + 10, Color0[mask_nm][i], (np.log10(StellarMass1[index[mask_nm]])[i] + 10)-(np.log10(StellarMass0[mask_nm])[i] + 10), Color1[index[mask_nm]][i]-Color0[mask_nm][i], ls='--', color='gray', length_includes_head=True, width=0.0005 )

plt.xlim( [9.0,12] )
plt.ylim( [0.1,0.8] )
plt.xlabel( '$\log M_{\star} [10^{10}\ M_{\odot}]$', fontsize = 18 )
plt.ylabel( 'g - r', fontsize = 18 )




plt.show()

