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

#========================================================================================
#   PARAMETERS
#========================================================================================
#Data folder
DataFolder = '/home/bustamsn/PhD/Data/cosmological_BH/'
#DataFolder = '/home/bustamsn/PhD/DataHaswell/'
Resolution = 256
#DataFolder = '/home/bustamsn/PhD/Spinstractor/data/'
DataResults = '.'
Simulations = ['cosmobh01', 'cosmobh02', 'cosmobh02']
Simulations = ['cosmobh03']
#Snapshot number
snapshot = 11
#Total bins
total_bins = 20
#Interpolation type
interptype = 'linear'


#Parameter of vonMisses function. (Dispersion)
K = [0,2,10]
colors = ['blue', 'green', 'red', 'black']
#========================================================================================
#   DATA
#========================================================================================

#Post processing

#On-fly data
Properties_onfly = np.loadtxt('%s/Sims%d/%s/analysis/BH_IDs.txt'%(DataFolder,Resolution,'cosmobh00'))
Properties_onfly = Properties_onfly[Properties_onfly[:,2] == 1.0]

#Concatenating data from raw files
#for k in K:
    #for i in Properties_onfly[:,0]:
        #os.system( 'tail -n 4 %s/Sims%3d/%s/analysis/analytical/spin_%d_k%d.txt >> %s/Sims%3d/%s/analysis/BH_IDs_%d.txt'%\
        #(DataFolder,Resolution,'cosmobh00',int(i),k,DataFolder,Resolution,'cosmobh00',k) )

plt.figure( figsize=(8,6) )
plt.subplots_adjust( top=0.9 )
ax1 = plt.subplot(1,1,1)
for i, Simulation in enumerate(Simulations):
    
    if i==1:
        None
        #for k in [K[i]]:
            #data = np.loadtxt('%s/Sims%3d/%s/analysis/BH_IDs_%d.txt'%(DataFolder,Resolution,'cosmobh00',k))
            #SpinBH = data[:,1]
            #MassBH = np.log10(data[:,2])+10
            
            ##SpinBH[MassBH>8] = np.random.random( sum(MassBH>8) )
            
            ##Median
            #bins, delta, median, perc25, perc75, std = Stats( MassBH, abs(SpinBH), total_bins = total_bins )
            #mask = (bins[1:]>6.5)&(bins[1:]<7)==False
            #median_f_bh1 = interp1d(bins[1:][mask], median[1:][mask], kind = interptype, bounds_error=False, fill_value='extrapolate')
            #perc25_f_bh1 = interp1d(bins[1:], perc25[1:], kind = interptype, bounds_error=False, fill_value='extrapolate')
            #perc75_f_bh1 = interp1d(bins[1:], perc75[1:], kind = interptype, bounds_error=False, fill_value='extrapolate')
            #bin_Mbh = np.log10(np.logspace( 5.9, 8.6, 1000 ))

            #mask_pos = median_f_bh1(bin_Mbh)>0
            #plt.plot(bin_Mbh[mask_pos], median_f_bh1(bin_Mbh)[mask_pos],'-',color=colors[i], lw=3, label = '$k=%d$'%K[i])
            #plt.fill_between( bin_Mbh[mask_pos], perc25_f_bh1(bin_Mbh)[mask_pos], perc75_f_bh1(bin_Mbh)[mask_pos], color=colors[i], alpha = 0.2 )
            #plt.plot(bin_Mbh[mask_pos], perc25_f_bh1(bin_Mbh)[mask_pos],color=colors[i], lw=1, ls = '-')
            #plt.plot(bin_Mbh[mask_pos], perc75_f_bh1(bin_Mbh)[mask_pos],color=colors[i], lw=1, ls = '-')
    
    else:
    
        SpinBH = np.array([])
        MassBH = np.array([])
        ModeBH = np.array([])
        
        for snapshot in [10,11]:
            SpinBH = np.append(SpinBH,np.loadtxt( '%s/Sims%3d/%s/analysis/BH_snap_%d/%s.txt'%(DataFolder,Resolution,Simulation,snapshot,'BH_SpinParameter') ))
            MassBH = np.append(MassBH,np.log10(np.loadtxt( '%s/Sims%3d/%s/analysis/BH_snap_%d/%s.txt'%(DataFolder,Resolution,Simulation,snapshot,'BH_Mass') ))+10)
            ModeBH = np.append(ModeBH,np.loadtxt( '%s/Sims%3d/%s/analysis/BH_snap_%d/%s.txt'%(DataFolder,Resolution,Simulation,snapshot,'BH_SpinModel') ))

        #Inf masks
        infmask = (np.isinf(MassBH) + np.isinf(SpinBH))==False
        MassBH = MassBH[infmask]
        SpinBH = SpinBH[infmask]
        ModeBH = ModeBH[infmask]
        mask_sort = np.argsort(MassBH)
        MassBH = MassBH[mask_sort]
        SpinBH = SpinBH[mask_sort]
        ModeBH = ModeBH[mask_sort]
        
        #Median
        plt.scatter( MassBH, abs(SpinBH), c = ModeBH, lw=0 )
        
        #bins, delta, median, perc25, perc75, std = Stats( MassBH, abs(SpinBH), total_bins = total_bins )
        #median_f_bh1 = interp1d(bins[1:], median[1:], kind = interptype, bounds_error=False, fill_value='extrapolate')
        #perc25_f_bh1 = interp1d(bins[1:], perc25[1:], kind = interptype, bounds_error=False, fill_value='extrapolate')
        #perc75_f_bh1 = interp1d(bins[1:], perc75[1:], kind = interptype, bounds_error=False, fill_value='extrapolate')
        #bin_Mbh = np.log10(np.logspace( 5.9, 8.6, 1000 ))

        #mask_pos = median_f_bh1(bin_Mbh)>0
        #plt.plot(bin_Mbh[mask_pos], median_f_bh1(bin_Mbh)[mask_pos],'-',color=colors[i], lw=3, label = '$k=%d$'%K[i])
        #plt.fill_between( bin_Mbh[mask_pos], perc25_f_bh1(bin_Mbh)[mask_pos], perc75_f_bh1(bin_Mbh)[mask_pos], color=colors[i], alpha = 0.2 )
        #plt.plot(bin_Mbh[mask_pos], perc25_f_bh1(bin_Mbh)[mask_pos],color=colors[i], lw=1, ls = '-')
        #plt.plot(bin_Mbh[mask_pos], perc75_f_bh1(bin_Mbh)[mask_pos],color=colors[i], lw=1, ls = '-')


datax, datay, datax0, datax1, datay0, datay1 = np.loadtxt( "./data/Reynolds_spin.txt" ).T
datax0 = abs(np.log10((datax-datax0)*1e6) - np.log10((datax)*1e6))
datax1 = abs(np.log10((datax-datax1)*1e6) - np.log10((datax)*1e6))
plt.errorbar(np.log10(datax*1e6), datay, xerr = [datax0, datax1], yerr = [datay0, datay1], 
             fmt='D', lw = 2, ecolor='black', color='orange', ms=8, label = 'Reynolds 2013' )

#Formatting plot
plt.ylim( (-0.05, 1.05) )
plt.xlim( (5.68, 8.7) )  
plt.legend( loc = 'lower central', ncol=4, bbox_to_anchor=(1.03, 1.12), fontsize = 15 )
plt.xlabel('$\log_{10}\ \mathrm{M}_{\mathrm{bh}}\ [\mathrm{M}_{\odot}]$', fontsize = 20)
plt.ylabel('$|a|$', fontsize = 20)
ax1.tick_params(labelsize=18)
#plt.savefig('./figs/SpinDistro.pdf')

plt.show()
