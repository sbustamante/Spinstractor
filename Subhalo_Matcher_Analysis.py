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
Simulation1 = 'cosmobh04'

#========================================================================================
#   PLOTS
#========================================================================================
StellarMass0 = np.loadtxt( './data/stellarmass_%s.txt'%Simulation0 )
StellarMass1 = np.loadtxt( './data/stellarmass_%s.txt'%Simulation1 )

#StellarHalfRadius0 = np.loadtxt( './data/stellarhalfradius_%s.txt'%Simulation0 )
#StellarHalfRadius1 = np.loadtxt( './data/stellarhalfradius_%s.txt'%Simulation1 )

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
residual2 = np.zeros( len(StellarMass0) )
ID1 = np.arange( len(StellarMass1) )
for i,r in enumerate(Spos0):
    index[i] = int(ID1[np.argsort( np.linalg.norm( Spos1 - r, axis=1 ) )[0]])
    residual[i] = np.sort( np.linalg.norm( Spos1 - r, axis=1 ) )[0]
    residual2[i] = np.sort( np.linalg.norm( Spos1 - r, axis=1 ) )[1]
    
index = index.astype(int)
    
#Discarding not-matching haloes
mask_nm = (residual<=1e2)#*(np.log10(StellarMass1[index])+10>9.5)#*(Color0>0.6)


plt.figure( figsize=(15,18) )
plt.subplots_adjust( top = 0.94, left = 0.06, bottom = 0.06, right = 0.98 )

plt.subplot( 4,2,1 )
plt.title( 'Fiducial' )
plt.scatter( np.log10(StellarMass0) + 10 , Color0, alpha = 0.5, marker = 'o', s=20, c=DT0, vmin=0, vmax=1, cmap = 'jet_r' )
plt.xlim( [8.5,12] )
plt.ylim( [0.0,1.0] )
cb1 = plt.colorbar()
plt.xlabel( '$\log M_{\star} [10^{10}\ M_{\odot}]$', fontsize = 18 )
plt.ylabel( 'g - r', fontsize = 18 )
cb1.set_label( 'D / T', fontsize = 18 )

plt.subplot( 4,2,2 )
plt.title( 'SG switch' )
plt.scatter( np.log10(StellarMass1) + 10 , Color1, alpha = 0.5, marker = 'o', s=20, c=DT1, vmin=0, vmax=1, cmap = 'jet_r' )
plt.xlim( [8.5,12] )
plt.ylim( [0.0,1.0] )
cb1 = plt.colorbar()
plt.xlabel( '$\log M_{\star} [10^{10}\ M_{\odot}]$', fontsize = 18 )
plt.ylabel( 'g - r', fontsize = 18 )
cb1.set_label( 'D / T', fontsize = 18 )


plt.subplot(4,2,3)
plt.plot( np.log10(StellarMass0[mask_nm])+10, np.log10(StellarMass1[index[mask_nm]])-np.log10(StellarMass0[mask_nm]), 'o', alpha=0.5 )
#plt.plot( np.log10(StellarMass0[mask_nm])+10, np.log10(DMMass1[index[mask_nm]])-np.log10(DMMass0[mask_nm]), 'o', alpha=0.5 )
plt.xlabel( '$\log(M_{\star, fid}/10^{10} M_{\odot})$', fontsize = 20 )
plt.ylabel( '$\Delta\logM_{\star}$', fontsize = 20 )
#plt.ylim( (-0.5,0.5) )
plt.xlim( (7,12) )
plt.hlines( 0, 7, 12, linestyles='--', color = 'grey', lw = 2 )


plt.subplot(4,2,4)
#plt.plot( np.log10(StellarMass0[mask_nm])+10, Color1[index[mask_nm]]-Color0[mask_nm], 'o', alpha=0.5 )



mask_nm = (residual<=1e2)*(np.log10(StellarMass1[index])+10>9.5)
#plt.plot(  Color0[mask_nm], DT0[mask_nm], 'o', alpha=0.5 )
#plt.plot(  Color1[index[mask_nm]], DT1[index[mask_nm]], 'o', alpha=0.5, color='red' )

Counts = np.histogram2d( Color0[mask_nm], DT0[mask_nm], range = ((0.1,0.8),(0,1)), bins = 6)[0]
plt.contour( np.transpose(Counts), extent = (0.1,0.8,0,1), antialiased=True, colors='blue')

Counts = np.histogram2d( Color1[index[mask_nm]], DT1[index[mask_nm]], range = ((0.1,0.8),(0,1)), bins = 6)[0]
plt.contour( np.transpose(Counts), extent = (0.1,0.8,0,1), antialiased=True, colors='red')
mask_nm = (residual<=1e2)




plt.xlabel( '$\log(M_{\star, fid}/10^{10} M_{\odot})$', fontsize = 20 )
plt.ylabel( '$\Delta\ \mathrm{g - r}$', fontsize = 20 )
#plt.xlim( (7,12) )
#plt.hlines( 0, 7, 12, linestyles='--', color = 'grey', lw = 2 )


plt.subplot(4,2,5)
plt.plot( np.log10(StellarMass0[mask_nm])+10, DT1[index[mask_nm]]-DT0[mask_nm], 'o', alpha=0.5 )
plt.xlabel( '$\log(M_{\star, fid}/10^{10} M_{\odot})$', fontsize = 20 )
plt.ylabel( '$\Delta\ \mathrm{D / T}$', fontsize = 20 )
plt.xlim( (7,12) )
plt.hlines( 0, 7, 12, linestyles='--', color = 'grey', lw = 2 )


plt.subplot(4,2,6)
#plt.plot( np.log10(StellarMass0[mask_nm])+10, GasMass1[index[mask_nm]]/DMMass1[index[mask_nm]]-GasMass0[mask_nm]/DMMass0[mask_nm], 'o', alpha=0.5 )
plt.xlabel( '$\log(M_{\star, fid}/10^{10} M_{\odot})$', fontsize = 20 )
plt.ylabel( '$\Delta\ M_{\mathrm{gas}}/M_{\mathrm{DM}}$', fontsize = 20 )
plt.xlim( (7,12) )
plt.hlines( 0, 7, 12, linestyles='--', color = 'grey', lw = 2 )


plt.subplot(4,2,7)
mask_nm = (residual<=1e2)*(np.log10(StellarMass1[index])+10>9.5)*(Color0<0.6)
#mask_nm = (residual<=1e2)*(np.log10(StellarMass1[index])+10>9.5)
#Histogram
Counts = np.histogram2d( Color1[index[mask_nm]]-Color0[mask_nm], DT1[index[mask_nm]]-DT0[mask_nm], range = ((-0.4,0.4),(-0.6,0.6)), bins = 4)[0]
plt.contour( np.transpose(Counts), extent = (-0.3,0.3,-0.6,0.6), antialiased=True, colors='black')#, levels = [2,6,10] )
plt.hlines( 0, -0.6, 0.6, linestyles='--', color = 'grey', lw = 2 )
plt.vlines( 0, -1, 1, linestyles='--', color = 'grey', lw = 2 )

#plt.plot( Color1[index[mask_nm]]-Color0[mask_nm], DT1[index[mask_nm]]-DT0[mask_nm], 'o', alpha=0.5 )
plt.scatter( Color1[index[mask_nm]]-Color0[mask_nm], DT1[index[mask_nm]]-DT0[mask_nm], c=np.log10(StellarMass1[index[mask_nm]])-np.log10(StellarMass0[mask_nm]), marker='o', alpha=0.5, s=30, cmap='seismic', vmin=-0.7, vmax=0.7 )
cb1 = plt.colorbar()
plt.xlabel( '$\Delta\ \mathrm{g - r}$', fontsize = 20 )
plt.ylabel( '$\Delta\ \mathrm{D / T}$', fontsize = 20 )
cb1.set_label( '$\Delta log M_{\star}$', fontsize = 18 )
plt.xlim( (-0.3,0.3) )
plt.ylim( (-0.6,0.6) )



#plt.subplot(4,2,8)
#mask_nm = (residual<=1e2)*(np.log10(StellarMass1[index])+10>9.5)#*(Color0>0.6)
##Histogram
#Counts = np.histogram2d( Color1[index[mask_nm]]-Color0[mask_nm], DT1[index[mask_nm]]-DT0[mask_nm], range = ((-0.4,0.4),(-0.6,0.6)), bins = 4)[0]
#plt.contour( np.transpose(Counts), extent = (-0.3,0.3,-0.6,0.6), antialiased=True, colors='black', levels = [2,6,10] )
#plt.hlines( 0, -0.6, 0.6, linestyles='--', color = 'grey', lw = 2 )
#plt.vlines( 0, -1, 1, linestyles='--', color = 'grey', lw = 2 )

##plt.plot( Color1[index[mask_nm]]-Color0[mask_nm], DT1[index[mask_nm]]-DT0[mask_nm], 'o', alpha=0.5 )
#plt.scatter( Color1[index[mask_nm]]-Color0[mask_nm], (StellarHalfRadius1[index[mask_nm]]-StellarHalfRadius0[mask_nm])/StellarHalfRadius0[mask_nm], c=DT1[index[mask_nm]]-DT0[mask_nm], marker='o', alpha=0.5, s=30, cmap='seismic')#, vmin=-0.7, vmax=0.7 )
#cb1 = plt.colorbar()
#plt.xlabel( '$\Delta\ \mathrm{g - r}$', fontsize = 20 )
#plt.ylabel( '$\Delta\ \mathrm{D / T}$', fontsize = 20 )
#cb1.set_label( '$\Delta log M_{\star}$', fontsize = 18 )
##plt.xlim( (-0.3,0.3) )
##plt.ylim( (-0.6,0.6) )



#plt.subplot( 4,2,8 )
#mask_nm = (residual<=1e2)*(np.log10(StellarMass1[index])+10>9.5)*(Color0>0.6)
#plt.scatter( np.log10(StellarMass0[mask_nm]) + 10 , Color0[mask_nm], alpha = 0.5, marker = 's', s=20, c=DT0[mask_nm], vmin=0, vmax=1, cmap = 'jet_r' )
#plt.scatter( np.log10(StellarMass1[index[mask_nm]]) + 10 , Color1[index[mask_nm]], alpha = 0.5, marker = 'o', s=20, c=DT1[index[mask_nm]], vmin=0, vmax=1, cmap = 'jet_r' )

#for i in xrange( sum(mask_nm) ):
    #plt.plot( [np.log10(StellarMass0[mask_nm])[i] + 10, np.log10(StellarMass1[index[mask_nm]])[i] + 10], [Color0[mask_nm][i], Color1[index[mask_nm]][i]], '--', color='gray' )

#plt.xlim( [9.9,12] )
#plt.ylim( [0.6,0.8] )
#cb1 = plt.colorbar()
#plt.xlabel( '$\log M_{\star} [10^{10}\ M_{\odot}]$', fontsize = 18 )
#plt.ylabel( 'g - r', fontsize = 18 )
#cb1.set_label( 'D / T', fontsize = 18 )




plt.show()

