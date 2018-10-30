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
from matplotlib import rc
rc('font',family='Palatino')
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#==================================================================================================
#       FUNCTIONS
#==================================================================================================
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
    #return bins[:-1], delta, median, perc25, perc75, std
    return np.insert(0.5*(bins[1:] + bins[:-1]), 0, bins[0]), delta, np.insert(median, 0, median[0]), np.insert(perc25, 0, perc25[0]), np.insert(perc75, 0, perc75[0]), np.insert(std, 0, std[0])


#==================================================================================================
#       INPUT DATA
#==================================================================================================

#Data folder
#DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/'
#DataFolder = '/home/bustamsn/PhD/DataHaswell/'
DataFolder = '/home/bustamsn/PhD/Data/cosmological_BH/'
#Resolution
Resolution = 256
#Simulation
Simulations = ['cosmobh00', 'cosmobh01', 'cosmobh02', 'cosmobh03', 'cosmobh04']

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

    for i, snapshot in enumerate(np.arange(11,16)):
        
        MassBH = []
        MassHost = []
        
        #Loading subfind halos
        sub = arepo.Subfind('%s/Sims%3d/%s/output/groups_015'%(DataFolder, Resolution, sim), snapshot=snapshot, combineFiles=True )
        
        MassHost.append( np.log10(sub.SubhaloMassInRadType[:,4]*1e10) )
        MassBH.append( np.log10(sub.SubhaloBHMass*1e10) )
            
        MassBH = np.array(list(chain.from_iterable(MassBH)))
        MassHost = np.array(list(chain.from_iterable(MassHost)))
        
        #Inf masks
        infmask = (np.isinf(MassBH) + np.isinf(MassHost))==False
        MassBH = MassBH[infmask]
        MassHost = MassHost[infmask]
        
        
        bins, delta, median, perc25, perc75, std = Stats( MassHost, MassBH, total_bins = total_bins )
        median_f = interp1d(bins, median, kind = interptype, bounds_error=False, fill_value='extrapolate')
        perc25_f = interp1d(bins, perc25, kind = interptype, bounds_error=False, fill_value='extrapolate')
        perc75_f = interp1d(bins, perc75, kind = interptype, bounds_error=False, fill_value='extrapolate')
        bin_Mhost = np.linspace( 8.5, 11.2, 100 )
        
        plt.subplot( 5, 2, 2*i_sim+1 )
        plt.title( sim )
        plt.plot(bin_Mhost, median_f(bin_Mhost),'-',color=colors[i], lw=1, alpha=.8)
        #plt.fill_between( bin_Mhost, perc25_f(bin_Mhost), perc75_f(bin_Mhost), color=colors[i], alpha = 0.2, label = "z = %1.2f"%(1/sub.Time - 1) )
        plt.plot( Mbulg, Mbh, '--', lw = 1, label='Kormendy + Ho (2013)', color='red' )
        plt.xlim( (8.5, 11.1) )
        plt.ylim( (5, 10) )
        
        plt.subplot( 5, 2, 2*i+2 )
        plt.title( "z = %1.2f"%(1/sub.Time - 1) )
        plt.plot(bin_Mhost, median_f(bin_Mhost),'-',color=colors[i_sim], lw=1, alpha=.8)
        #plt.fill_between( bin_Mhost, perc25_f(bin_Mhost), perc75_f(bin_Mhost), color=colors[i_sim], alpha = 0.2, label = sim )
        plt.plot( Mbulg, Mbh, '--', lw = 1, label='Kormendy + Ho (2013)', color='red' )
        plt.xlim( (8.5, 11.1) )
        plt.ylim( (5, 10) )
        
        #plt.scatter( MassHost, MassBH, color='black', marker='.')

        
    plt.plot( Mbulg, Mbh, '--', lw = 1, label='Kormendy + Ho (2013)', color='red' )
    plt.xlim( (8.5, 11.1) )
    plt.ylim( (5, 10) )
    #plt.legend( loc='upper left', fancybox=True )
    plt.xlabel( '$\log_{10}(M_{\star}/M_{\odot} )$', fontsize = 12 )
    plt.ylabel( '$\log_{10}(M_{BH}/M_{\odot} )$', fontsize = 12 )
    #plt.legend( loc = 'upper left', fontsize = 8 )


#plt.savefig("figs/BHMass-HostMass_relation.pdf")
    
plt.show()
 
