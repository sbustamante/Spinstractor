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

#========================================================================================
#   FUNCTIONS
#========================================================================================
def Stats( X, Y, total_bins = 7 ):
    Y = Y[np.argsort(X)]
    X = X[np.argsort(X)]
    bins = np.array([ X[int(len(X)/total_bins*i)] for i in xrange(total_bins) ] + [X.max()])
    bins = np.linspace( X.min(), X.max(), total_bins+1 )
    delta = bins[1]-bins[0]
    idx  = np.digitize(X,bins)
    median = -1*np.ones( total_bins )
    std = -1*np.ones( total_bins )
    perc25 = -1*np.ones( total_bins )
    perc75 = -1*np.ones( total_bins )
    for k in range(1,total_bins+1):
        try:
            perc25[k-1] = np.percentile(Y[idx==k], 16)
            perc75[k-1] = np.percentile(Y[idx==k], 84)
            median[k-1] =  np.median(Y[idx==k])
            std[k-1] =  Y[idx==k].std()
        except:
            None
    return np.insert(0.5*(bins[1:] + bins[:-1]), 0, bins[0]), delta, np.insert(median, 0, median[0]), np.insert(perc25, 0, perc25[0]), np.insert(perc75, 0, perc75[0]), np.insert(std, 0, std[0])

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


#========================================================================================
#   PARAMETERS
#========================================================================================
#Data folder
#DataFolder = '/home/bustamsn/PhD/Data/cosmological_BH/'
#DataFolder = '/home/bustamsn/PhD/DataHaswell/'
DataFolder = '/home/bustamsn/PhD/DataFreya/'
Resolution = 512
#DataFolder = '/home/bustamsn/PhD/Spinstractor/data/'
DataResults = '.'
#Simulations = ['cosmobh00', 'cosmobh02', 'cosmobh03']
Simulations = ['cosmobh03']
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
for i, Simulation in enumerate(Simulations):
    
        dMassBH = np.array([])
        SpinBH = np.array([])
        dMassBH_Edd = np.array([])
        MassBH = np.array([])
        ModeBH = np.array([])
        ModeRM = np.array([])
        
        for j, snapshot in enumerate([6, 10, 12, 15]):
            SpinBH = np.append(SpinBH,np.loadtxt( '%s/Sims%3d/%s/analysis/BH_snap_%d/%s.txt'%(DataFolder,Resolution,Simulation,snapshot,'BH_SpinParameter') ))
            dMassBH = np.append(dMassBH,np.loadtxt( '%s/Sims%3d/%s/analysis/BH_snap_%d/%s.txt'%(DataFolder,Resolution,Simulation,snapshot,'BH_MdotBondi') ))
            dMassBH_Edd = np.append(dMassBH_Edd,np.loadtxt( '%s/Sims%3d/%s/analysis/BH_snap_%d/%s.txt'%(DataFolder,Resolution,Simulation,snapshot,'BH_MdotEddington') ))
            MassBH = np.append(MassBH,np.log10(np.loadtxt( '%s/Sims%3d/%s/analysis/BH_snap_%d/%s.txt'%(DataFolder,Resolution,Simulation,snapshot,'BH_Mass') ))+10)
            ModeBH = np.append(ModeBH,np.loadtxt( '%s/Sims%3d/%s/analysis/BH_snap_%d/%s.txt'%(DataFolder,Resolution,Simulation,snapshot,'BH_SpinModel') ))
            ModeRM = np.append(ModeRM,np.loadtxt( '%s/Sims%3d/%s/analysis/BH_snap_%d/%s.txt'%(DataFolder,Resolution,Simulation,snapshot,'BH_CumMassGrowth_RM') ))>0

            #Sort masks
            mask_sort = np.argsort(MassBH)
            MassBH = MassBH[mask_sort]
            dMassBH = dMassBH[mask_sort]
            dMassBH_Edd = dMassBH_Edd[mask_sort]
            SpinBH = SpinBH[mask_sort]
            ModeBH = ModeBH[mask_sort]
            ModeRM = ModeRM[mask_sort].astype(int)
            
            Mdot = np.log10(dMassBH/dMassBH_Edd)
            if i>0:
                Mdot = np.log10(dMassBH/m_dot_eddignton(10**(MassBH-10), SpinBH))
            
            #Median
            plt.scatter( MassBH, Mdot, c = ModeBH, lw = 0, s = 50, alpha = 0.3 )
            #plt.scatter( MassBH, Mdot, c = colors[i], lw = 0, s = 50 )
            bins, delta, median, perc25, perc75, std = Stats( MassBH, Mdot, total_bins = total_bins[j] )
            median_f_bh1 = interp1d(bins[1:], median[1:], kind = interptype, bounds_error=False, fill_value='extrapolate')
            perc25_f_bh1 = interp1d(bins[1:], perc25[1:], kind = interptype, bounds_error=False, fill_value='extrapolate')
            perc75_f_bh1 = interp1d(bins[1:], perc75[1:], kind = interptype, bounds_error=False, fill_value='extrapolate')
            bin_Mbh = np.log10(np.logspace( 5.9, 8.6, 1000 ))

            #plt.plot(bin_Mbh, median_f_bh1(bin_Mbh),'-',color=colors[i], lw=3)
            #plt.fill_between( bin_Mbh, perc25_f_bh1(bin_Mbh), perc75_f_bh1(bin_Mbh), color=colors[i], alpha = 0.2 )
            #plt.plot(bin_Mbh, perc25_f_bh1(bin_Mbh),color=colors[i], lw=1, ls = '-')
            #plt.plot(bin_Mbh, perc75_f_bh1(bin_Mbh),color=colors[i], lw=1, ls = '-')


#Formatting plot
Mbharray = np.linspace( 5.5,10.0,200 )
Mbh_dot = np.linspace( -7, 0, 200 )

X0 = 0.002
beta = 2.0
Thr_feedback = np.log10(np.array([ np.min( (X0*(10**(m-8))**beta, 0.1) ) for m in Mbharray ]))

Mbh_thr1 = Mass_threshold_D( 10**Mbh_dot, 0.998 )
Mbh_thr2 = Mass_threshold_D( 10**Mbh_dot, 0.01 )


plt.plot( Mbharray, Thr_feedback, '--', color='black', lw = 3 )
plt.plot( np.log10(Mbh_thr1), Mbh_dot, '--', color='red', lw = 3 )
#plt.plot( np.log10(Mbh_thr2), Mbh_dot, '--', color='red', lw = 3 )
plt.hlines( -3, 5.5, 10, linestyles=':', color='black', lw = 3 )

plt.xlim( (5.5, 10.0) )
plt.ylim( (-7, 0) )
plt.legend( loc = 'lower central', ncol=4, bbox_to_anchor=(1.03, 1.12), fontsize = 15 )
plt.xlabel('$\log_{10}\ \mathrm{M}_{\mathrm{bh}}\ [\mathrm{M}_{\star}]$', fontsize = 20)
plt.ylabel('$\dot\mathrm{M}_{\mathrm{Bondi}}/\dot\mathrm{M}_{\mathrm{Edd}}$', fontsize = 20)
ax1.tick_params(labelsize=18)
plt.savefig('./figs/EvolutionDiagram1.pdf')

plt.show()
