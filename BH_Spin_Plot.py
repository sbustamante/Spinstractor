#========================================================================================
#   LIBRARIES
#========================================================================================
import matplotlib
matplotlib.use('Agg')
import numpy, h5py, IPython,astropy, pylab
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
import numpy as np
import matplotlib.pyplot as plt
#import BH_Physics as BHP
execfile( 'BH_Physics.py' )
import sys
plt.close('all')

#========================================================================================
#   PARAMETERS
#========================================================================================
#Data folder
DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/Sims256/'
#Simulation
Simulation = 'cosmobh00'
#Number of chunks (same number of used processors)
N_proc = 256
DataResults = '.'


#========================================================================================
#   LIBRARIES
#========================================================================================
BH = black_hole_sim( simulation = Simulation, snapbase = 'snapshot', n_snap = 1000, datafolder = DataFolder, resultsfolder = DataResults )
indexes = np.loadtxt('%s%s/analysis/BH_IDs.txt'%(DataFolder,Simulation))

for i in indexes[indexes[:,2]==1,0].astype(int):

    #Calculating merger tree correction
    BH.merger_tree_correction(i)
    ndata = BH.loading_spin( i )

    plt.figure( figsize = (10,9) )
    plt.subplots_adjust( hspace = 0.1 )
    plt.rc('font', size=18)


    Time_inter = np.linspace( BH.t.min(), BH.t.max(), 1000 )
    time = Time_inter
    costh = np.sum( BH.Lgas(Time_inter)*BH.SpinDir(Time_inter), axis = 0 )

    #Spin evolution
    plt.subplot(4,1,1)
    plt.plot( time, abs(BH.A(Time_inter)), '-', zorder = 10, lw = 1.5, color='black' )
    plt.pcolor( BH.t, [0,1.1], np.array([BH.mode]), cmap = 'jet', snap = True, alpha = 0.3 )
    plt.vlines( BH.Tmerger, 0, 1.02, linestyle='--', color='gray' )

    plt.xlim((time.min(), time.max()))
    plt.ylim([0,1.02])
    plt.ylabel( '$|a|$', fontsize = 26 )
    plt.legend( fontsize = 18, ncol=6, bbox_to_anchor=(0.95, 1.3) )
    plt.xticks([])


    #Radiative efficiency
    plt.subplot(4,1,2)
    plt.plot( time, BH.radiative_efficiency((BH.A(Time_inter))), '-', zorder = 10, lw = 1.5, color='black' )
    plt.pcolor( BH.t, [-0.05,0.4], np.array([BH.mode]), cmap = 'jet', snap = True, alpha = 0.3 )
    plt.ylim([-0.05,0.4])
    plt.xlim((time.min(), time.max()))
    plt.ylabel( '$\epsilon_r$', fontsize = 26 )
    plt.xticks([])
    plt.yticks( [0.0,0.1,0.2,0.3] )
    plt.vlines( BH.Tmerger, -0.05, 0.4, linestyle='--', color='gray' )

    #Spin orientation
    plt.subplot(4,1,3)
    plt.plot( time, costh, '-', zorder = 10, lw = 1.5, color='black' )
    plt.pcolor( BH.t, [-1,1.19], np.array([BH.mode]), cmap = 'jet', snap = True, alpha = 0.3 )
    plt.xlim((time.min(), time.max()))
    plt.ylim([-1,1.19])
    plt.xticks([])
    plt.ylabel( '$\cos(\\theta)$', fontsize = 26 )
    plt.vlines( BH.Tmerger, -1, 1.19, linestyle='--', color='gray' )

    #Accretion rate
    plt.subplot(4,1,4)
    plt.plot( time, np.log10(BH.Mbh(Time_inter))+10, '-', zorder = 10, color = 'black', lw = 1.5 )
    plt.pcolor( BH.t, [4,10], np.array([BH.mode]), cmap = 'jet', snap = True, alpha = 0.3 )
    plt.xlim((time.min(), time.max()))
    plt.ylim( (4,10) )
    plt.ylabel( 'M$_{bh}$ [M$_{\odot}$]' )
    plt.xlabel( 'time since formation [Gyr]', fontsize = 18 )
    plt.vlines( BH.Tmerger, 4, 10, linestyle='--', color='gray' )

        
        
    plt.savefig('./fig_%d.png'%(i))
