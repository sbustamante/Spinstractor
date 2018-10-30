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

#Observational Kormendi-Ho function
Mbulg = np.linspace( 8, 12, 100 )
Mbh = 1.06352*( Mbulg - np.log10(2.7e8) ) + 6

#========================================================================================
#   PARAMETERS
#========================================================================================
#Data folder
#DataFolder = '/home/bustamsn/PhD/Data/cosmological_BH/'
DataFolder0 = '/home/bustamsn/PhD/DataHaswell/'
DataFolder1 = '/home/bustamsn/PhD/DataFreya/'
Resolution = 512
#DataFolder = '/home/bustamsn/PhD/Spinstractor/data/'
DataResults = '.'
#Simulation
Sim0 = 'cosmobh00'
Sim1 = 'cosmobh03'
#Snapshot number
snapshot = 15
#Total bins
total_bins = 20
#Interpolation type
interptype = 'linear'


#Parameter of vonMisses function. (Dispersion)
K = [0,10,2]
colors = ['blue', 'red', 'green', 'black']
#========================================================================================
#   DATA
#========================================================================================

#Loading subfind subhalos
sub0_h = arepo.Subfind('%s/Sims%3d/%s/output/groups_015'%(DataFolder0, Resolution, Sim0), snapshot=snapshot, combineFiles=True )
sub1_h = arepo.Subfind('%s/Sims%3d/%s/output/groups_015'%(DataFolder1, Resolution, Sim1), snapshot=snapshot, combineFiles=True )

#Total bins
total_bins = 20
#Interpolation type
interptype = 'linear'

#Observational Kormendi-Ho function
Mbulg = np.linspace( 8, 12, 100 )
Mbh = 1.06352*( Mbulg - np.log10(2.7e8) ) + 6

#Fiducial model-------------------------------------------------
MassHost = np.log10(sub0_h.SubhaloMassInRadType[:,4]*1e10)
MassBH = np.log10(sub0_h.SubhaloBHMass*1e10)
#Inf masks
infmask = (np.isinf(MassBH) + np.isinf(MassHost))==False
MassBH = MassBH[infmask]
MassHost = MassHost[infmask]
#Median
bins, delta, median, perc25, perc75, std = Stats( MassHost, MassBH, total_bins = total_bins )
median_f_h = interp1d(bins, median, kind = interptype, bounds_error=False, fill_value='extrapolate')
perc25_f_h = interp1d(bins, perc25, kind = interptype, bounds_error=False, fill_value='extrapolate')
perc75_f_h = interp1d(bins, perc75, kind = interptype, bounds_error=False, fill_value='extrapolate')
bin_Mhost = np.linspace( 7, 11.2, 100 )

#Sub1 model-------------------------------------------------
MassHost = np.log10(sub1_h.SubhaloMassInRadType[:,4]*1e10)
MassBH = np.log10(sub1_h.SubhaloBHMass*1e10)
#Inf masks
infmask = (np.isinf(MassBH) + np.isinf(MassHost))==False
MassBH = MassBH[infmask]
MassHost = MassHost[infmask]
#Median
bins, delta, median, perc25, perc75, std = Stats( MassHost, MassBH, total_bins = total_bins )
median_f1_h = interp1d(bins, median, kind = interptype, bounds_error=False, fill_value='extrapolate')
perc25_f1_h = interp1d(bins, perc25, kind = interptype, bounds_error=False, fill_value='extrapolate')
perc75_f1_h = interp1d(bins, perc75, kind = interptype, bounds_error=False, fill_value='extrapolate')


#Formatting figure

plt.figure( figsize=(8,6) )
plt.subplots_adjust( top=0.9 )
ax1 = plt.subplot(1,1,1)


#plt.scatter( np.log10(sub0.SubhaloMassInRadType[:,4]*1e10), np.log10(sub0.SubhaloBHMass*1e10), color='red', marker='.', s=3)
#plt.scatter( np.log10(sub1.SubhaloMassInRadType[:,4]*1e10), np.log10(sub1.SubhaloBHMass*1e10), color='blue', marker='.', s=3)

plt.plot(bin_Mhost, median_f_h(bin_Mhost),'-',color='gray', lw=3, label = 'Fiducial')
#plt.fill_between( bin_Mhost, perc25_f_h(bin_Mhost), perc75_f_h(bin_Mhost), color='gray', alpha = 0.2 )
plt.plot(bin_Mhost, perc25_f_h(bin_Mhost),color='gray', lw=1, ls = '-')
plt.plot(bin_Mhost, perc75_f_h(bin_Mhost),color='gray', lw=1, ls = '-')

plt.plot(bin_Mhost, median_f1_h(bin_Mhost),'-',color='blue', lw=3, label = 'SG, $k=10$')
#plt.fill_between( bin_Mhost, perc25_f1_h(bin_Mhost), perc75_f1_h(bin_Mhost), color='blue', alpha = 0.2 )
plt.plot(bin_Mhost, perc25_f1_h(bin_Mhost),color='blue', lw=1, ls = '-')
plt.plot(bin_Mhost, perc75_f1_h(bin_Mhost),color='blue', lw=1, ls = '-')

#plt.scatter( np.log10(sub2.SubhaloMassInRadType[:,4]*1e10), np.log10(sub2.SubhaloBHMass*1e10), color='red', marker='.', s=3)

#plt.plot(bin_Mhost, median_f2_h(bin_Mhost),'-',color='red', lw=3, label = '$k=10$')
#plt.fill_between( bin_Mhost, perc25_f2_h(bin_Mhost), perc75_f2_h(bin_Mhost), color='red', alpha = 0.2 )
#plt.plot(bin_Mhost, perc25_f2_h(bin_Mhost),color='red', lw=1, ls = '-')
#plt.plot(bin_Mhost, perc75_f2_h(bin_Mhost),color='red', lw=1, ls = '-')

plt.plot( Mbulg, Mbh, '--', color='black', lw=4, label='Kormendi-Ho', alpha = 0.8 )
plt.xlim( (8, 11) )
plt.ylim( (5.8, 9.3) )
plt.legend( loc = 'lower central', ncol=4, bbox_to_anchor=(1.04, 1.12), fontsize = 15 )
plt.xlabel( '$\log_{10}\ \mathrm{M}_{\star}\ [\mathrm{M}_{\odot}]$', fontsize = 20 )
plt.ylabel( '$\log_{10}\ \mathrm{M}_{\mathrm{bh}}\ [\mathrm{M}_{\odot}]$', fontsize = 20 )
plt.xticks( fontsize=16 )
plt.yticks( fontsize=16 )
ax1.tick_params(labelsize=18)
plt.savefig('./figs/MassDistro2.pdf')

plt.show()
