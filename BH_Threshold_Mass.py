#========================================================================================
#   LIBRARIES
#========================================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import BH_Physics as BHP
import arepo
import os
import matplotlib
from matplotlib import rc
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.close('all')

#========================================================================================
#		GLOBAL VARIABLES
#========================================================================================
#Cavendish constant [SI]
GC	=	6.67384e-11
#Light speed
CLIGHT	=	2.99792458e8
#Proton Mass
PROTONMASS =    1.67262178e-27
#Thomson constant
THOMPSON =      6.65245873e-29
#1e8 Solar units
E8MSUN = 	1.989e38
#Seconds per year
SECPERYEAR =    3.15569e7
#Solar mass
SOLARMASS =     1.989e30


#Standard units
units = {'M':1.989e40, 'L':3.08568e19, 'T':3.15569e16, 'V':1000. }

def radiative_efficiency( a ):
    try:
        len(a)
        a[a > 0.998] = 0.998
    except:
        None
    
    Z1 = 1 + (1-a**2)**(1/3.0)*( (1+a)**(1/3.0) + (1-a)**(1/3.0) )
    Z2 = ( 3*a**2 + Z1**2 )**(0.5)
    Rlso = 3 + Z2 - abs(a)/a*( (3-Z1)*(3+Z1+2*Z2) )**(0.5)
    return 1 - np.sqrt( 1 - 2/(3*Rlso) )

def m_dot_eddignton( Mbh, Spin ):
    rad_eff = radiative_efficiency( Spin )
    return (4*np.pi*GC*CLIGHT*PROTONMASS/(rad_eff*CLIGHT**2*THOMPSON))*Mbh*units['T']

def Mass_threshold_D( Mdot, a ):
    #Parameters
    alpha = 0.1
    v2v1 = 2*( 1 + 7*alpha)/( alpha**2*(4 + alpha**2) )
    return 0.447228*abs(a)**(-27/47.)*radiative_efficiency(a)**(64/235.)*Mdot**(-2/47.)*(v2v1)**(27/47.)*alpha**(44/47.)*1e8

def Mass_threshold_W( Mdot, a ):
    #Parameters
    alpha = 0.1
    v2v1 = 2*( 1 + 7*alpha)/( alpha**2*(4 + alpha**2) )
    return ( 3.6e3*abs(a)**(5/8.)*Mdot**(-1/4.)*v2v1**(-5/8.)*alpha**(-1/2.) )**-8*1e8

def Mass_threshold_S( Mdot, a ):
    #Parameters
    alpha = 0.1
    v2v1 = 2*( 1 + 7*alpha)/( alpha**2*(4 + alpha**2) )
    return ( 1.5e3*radiative_efficiency(a)**(8/27.)*Mdot**(-8/27.)*alpha**(14/27.) )**(27/26)*1e8


#========================================================================================
#   PARAMETERS
#========================================================================================
#Data folder
DataFolder = '/home/bustamsn/PhD/Data/cosmological_BH/'
#DataFolder = '/home/bustamsn/PhD/DataHaswell/'
Resolution = 256
#DataFolder = '/home/bustamsn/PhD/Spinstractor/data/'
DataResults = '.'
Simulations = ['cosmobh00', 'cosmobh02', 'cosmobh03']
#Snapshot number
snapshot = 15
#Total bins
total_bins = [10,10,10,10]
#Interpolation type
interptype = 'linear'


#Parameter of vonMisses function. (Dispersion)
K = [0,10]
colors = ['gray', 'blue', 'red']
#========================================================================================
#   DATA
#========================================================================================

plt.figure( figsize=(8,6) )
plt.subplots_adjust( top=0.9 )
ax1 = plt.subplot(1,1,1)
#Formatting plot
Mbharray = np.linspace( 5.5,10.0,200 )
Mbh_dot = np.linspace( -6, 0, 200 )

X0 = 0.002
beta = 2.0
Thr_feedback = np.log10(np.array([ np.min( (X0*(10**(m-8))**beta, 0.1) ) for m in Mbharray ]))

for a in [ 0.001, 0.01, 0.1, 0.998 ]:
    Mbh_thr1 = Mass_threshold_D( 10**Mbh_dot, a )
    Mbh_thr2 = Mass_threshold_S( 10**Mbh_dot, a )
    plt.plot( np.log10(Mbh_thr1), Mbh_dot, '--', lw = 2,label = "$a=%1.3f$"%(a) )
    #plt.plot( np.log10(Mbh_thr2), Mbh_dot, '-', lw = 2,label = "$a=%1.3f$"%(a) )

plt.xlim( (5.5, 10.0) )
plt.ylim( (-6, 0) )
plt.legend( loc = 'lower central', ncol=1, fontsize = 15 )
plt.xlabel('$\log_{10}\ \mathrm{M}_{\mathrm{bh}}\ [\mathrm{M}_{\odot}]$', fontsize = 20)
plt.ylabel('$\dot\mathrm{M}_{\mathrm{Bondi}}/\dot\mathrm{M}_{\mathrm{Edd}}$', fontsize = 20)
ax1.tick_params(labelsize=18)
#plt.savefig('./figs/EvolutionDiagram2.png')

plt.show()
