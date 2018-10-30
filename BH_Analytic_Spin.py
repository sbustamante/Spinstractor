#========================================================================================
#   LIBRARIES
#========================================================================================
import numpy as np
import matplotlib.pyplot as plt
import BH_Physics as BHP
from scipy.interpolate import interp1d
from astropy.cosmology import WMAP9
import matplotlib
from matplotlib import rc
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
import sys
import os
plt.close('all')


#========================================================================================
#   DATA
#========================================================================================
#Data folder
DataFolder = '/home/bustamsn/PhD/Data/cosmological_BH/Sims256/'
#DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/'
#DataFolder = '/home/bustamsn/PhD/DataHaswell/Sims512/'
DataResults = '.'
Simulation = 'cosmobh00'
#Values of anisotropy to compute
K_values = [0,2,10]

#========================================================================================
#   SPIN HISTORIES
#========================================================================================
BH = BHP.black_hole_sim( simulation = Simulation, snapbase = 'snapshot', n_snap = 1000, datafolder = DataFolder, resultsfolder = DataResults )
indexes = np.loadtxt('%s%s/analysis/BH_IDs.txt'%(DataFolder,Simulation))
#Creating folder to store analytical profiles
os.system('mkdir %s%s/analysis/analytical'%(DataFolder,Simulation))

mask_z0 = indexes[:,2] == 1
N_files = int(np.sum(mask_z0))
#Argument 1 in console: number of current chunk
chunk = int(sys.argv[1])
#Argument 2 in console: total number of chunks in which data is divided (parallel computing)
total_chunks = int(sys.argv[2])

N0 = int(N_files*(chunk-1)/total_chunks)
N1 = int(N_files*(chunk)/total_chunks)
if chunk == total_chunks:
    N1 = N_files

for i in [518]:#indexes[mask_z0][N0:N1,0].astype(int):
    plt.close('all')
    print chunk, i
    BH.merger_tree_correction(i)
    ndata = BH.loading_spin( i )
    
    time = np.linspace( BH.t.min(), BH.t.max(), 10000 )
    redshift_range = np.insert(np.logspace( -3, 4, 10000), 0, 0)
    redshift_f = interp1d(WMAP9.lookback_time( redshift_range ), redshift_range)
    redshift = redshift_f(time[-1]-time)
    
    plt.figure(figsize = (12,10))
    plt.subplots_adjust( hspace=0.2, top = 0.94, right = 0.99, left = 0.08, bottom = 0.06 )
    #plt.subplot(5,1,1)
    #plt.plot(BH.t, abs(BH.a), lw = 2, color='black') 
    #plt.subplot(5,1,2)
    #plt.plot(BH.t, np.log10(BH.Mbh(BH.t)), lw = 2, color='black')
    #plt.subplot(5,1,3)
    #plt.plot(BH.t, BH.radiative_efficiency(BH.a), lw = 2, color='black')
    #plt.subplot(5,1,5)
    #plt.plot(BH.t, np.log10(BH.Mbh_dot(BH.t)/BH.m_dot_eddignton( BH.Mbh(BH.t), 0.1 )), lw = 2, color='black') 
    
    for k in K_values:
        if ndata>10:
            print "k=",k
            BH.merger_tree_correction(i)
            ndata = BH.loading_spin( i )
            BH.spin_evolution( mode='auto', rad_eff = -1, alpha = 0.1, kparam = k, a0 = 0.05, Nmax = 5e5 )
            
            #Plotting------------------------------------------------------------------------------
            #Coherent/Self-gravity regime shades
            if k == 0:
                shade_regime = interp1d( BH.tms, BH.modeps, bounds_error = False, kind='nearest', fill_value='extrapolate' )
            
            #Spin evolution
            ax1 = plt.subplot(5,1,1)
            fs = interp1d( BH.tms, abs(BH.aps), bounds_error = False, kind='nearest', fill_value='extrapolate' )
            plt.plot( redshift, fs(time) )
            
            #Black hole mass
            ax2 = plt.subplot(5,1,2)
            fs = interp1d( BH.tms, np.log10(BH.mbhs)+10, bounds_error = False, kind='nearest', fill_value='extrapolate' )
            plt.plot( redshift, fs(time), label='$k=%d$'%(k) )
            
            #Radiative efficiency
            ax3 = plt.subplot(5,1,3)
            fs = interp1d( BH.tms, BH.radiative_efficiency(BH.aps), bounds_error = False, kind='nearest', fill_value='extrapolate' )
            plt.plot( redshift, fs(time) )
            
            #Angle
            ax4 = plt.subplot(5,1,4)
            fs = interp1d( BH.tms, BH.costh, bounds_error = False, kind='nearest', fill_value='extrapolate' )
            plt.plot( redshift, fs(time) )
            
            #Lambda
            ax5 = plt.subplot(5,1,5)
            fs = interp1d( BH.tms, BH.mdot, bounds_error = False, kind='nearest', fill_value='extrapolate' )
            plt.plot( redshift, np.log10(fs(time)) )
            
            np.savetxt( '%s%s/analysis/analytical/spin_%d_k%d.txt'%(DataFolder,Simulation,i,k), np.array( [BH.tms, BH.aps, BH.Mbh(BH.tms), BH.modeps] ).T, fmt='%e %e %e %d' )
          
    #Formating
    #Spin evolution
    ax1 = plt.subplot(5,1,1)
    plt.pcolor( redshift, [-0.05,1.05], np.array([shade_regime(time)]), cmap = 'jet', snap = True, alpha = 0.2, rasterized=True )
    plt.ylim( -0.05, 1.05 )
    plt.xlim( (0, redshift_f(time[-1]-time[0])) )
    ax1.set_xticklabels( [], fontsize=16 )
    ax1.set_yticks( np.arange( 0,1.1,0.2 ) )
    ax1.set_yticklabels( np.arange( 0,1.1,0.2 ), fontsize=16 )
    ax1.set_ylabel( "$|a|$", fontsize = 18 )
    ax1.yaxis.set_label_coords(-0.05,0.5)
    
    ax1_in = ax1.twiny()
    ax1_in.set_xticks( np.arange( 0,11,1 ) )
    ax1_in.set_xticklabels( ["%1.1f"%lab for lab in WMAP9.lookback_time(np.arange( 0,11,1 )).value], fontsize = 16 )
    ax1_in.set_xlabel( "Lookback time [Gyr]", fontsize = 18 )
    ax1_in.set_xlim( (0, redshift_f(time[-1]-time[0])) )
    
    
    #Black hole mass
    ax2 = plt.subplot(5,1,2)
    plt.pcolor( redshift, [5,10], np.array([shade_regime(time)]), cmap = 'jet', snap = True, alpha = 0.2, rasterized=True )
    plt.ylim( 5, 10 )
    plt.xlim( (0, redshift_f(time[-1]-time[0])) )
    plt.legend( loc='upper right', fontsize=16, fancybox=True, ncol=3 )
    ax2.set_xticklabels( [], fontsize=16 )
    ax2.set_yticks( np.arange( 5,11,1 ) )
    ax2.set_yticklabels( np.arange( 5,11,1 ), fontsize=16 )
    ax2.set_ylabel( "$\log\ \mathrm{M}_{\mathrm{bh}}\ \ [\mathrm{M}_{\odot}]$", fontsize = 18 )
    ax2.yaxis.set_label_coords(-0.05,0.5)
    
    #Radiative efficiency
    ax3 = plt.subplot(5,1,3)
    plt.pcolor( redshift, [0,0.4], np.array([shade_regime(time)]), cmap = 'jet', snap = True, alpha = 0.2, rasterized=True )
    plt.ylim( 0, 0.4 )
    plt.xlim( (0, redshift_f(time[-1]-time[0])) )
    ax3.set_xticklabels( [], fontsize=16 )
    ax3.set_yticks( [0,0.1,0.2,0.3,0.4] )
    ax3.set_yticklabels( ['$0\%$', '$10\%$', '$20\%$', '$30\%$', '$40\%$'], fontsize=16 )
    ax3.set_ylabel( "$\epsilon_{r}$", fontsize = 18 )
    ax3.yaxis.set_label_coords(-0.05,0.5)
    
    #Angle
    ax4 = plt.subplot(5,1,4)
    plt.pcolor( redshift, [-1.05,1.05], np.array([shade_regime(time)]), cmap = 'jet', snap = True, alpha = 0.2, rasterized=True )
    plt.ylim( -1.05, 1.05 )
    plt.xlim( (0, redshift_f(time[-1]-time[0])) )
    ax4.set_xticklabels( [], fontsize=16 )
    ax4.set_yticks( np.arange( -1,1.1,0.5 ) )
    ax4.set_yticklabels( np.arange( -1,1.1,0.5 ), fontsize=16 )
    ax4.set_ylabel( "$\cos\ \\theta$", fontsize = 18 )
    ax4.yaxis.set_label_coords(-0.05,0.5)
    
    #Lambda
    ax5 = plt.subplot(5,1,5)
    plt.pcolor( redshift, [-6,0.05], np.array([shade_regime(time)]), cmap = 'jet', snap = True, alpha = 0.2, rasterized=True )
    plt.ylim( -6, 0.05 )
    plt.hlines( -3, 0, redshift_f(time[-1]-time[0]), linestyles='--', color='gray'  )
    ax5.set_yticks( np.arange( -6,1,1 ) )
    ax5.set_yticklabels( np.arange( -6,1,1 ), fontsize=16 )
    ax5.set_ylabel( "$\log\ \chi$", fontsize = 18 )
    ax5.yaxis.set_label_coords(-0.05,0.5)
    ax5.set_xticks( np.arange( 0,11,1 ) )
    ax5.set_xticklabels( np.arange( 0,11,1 ), fontsize=16 )
    ax5.set_xlabel( "z", fontsize = 18 )
    ax5.set_xlim( (0, redshift_f(time[-1]-time[0])) )
          
    #plt.savefig('/home/bustamsn/PhD/CosmologicalBlackHoles/Spinstractor/figs/histories/BH_%d.png'%i)
    plt.savefig('./BH_%d.pdf'%i)
plt.show()
