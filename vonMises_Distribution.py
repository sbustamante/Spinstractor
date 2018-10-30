import numpy as np
import scipy as sp
import pylab as plt
from scipy.interpolate import interp1d
import os
import matplotlib
from matplotlib import rc
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

#from matplotlib import rc
#rc('font',family='Palatino')
#rc('font',**{'family':'serif','serif':['Palatino']})
#plt.close('all')
#mpl.rcParams['axes.linewidth'] = 2 #set the value globally



def von_Mises_PDF( CosTheta, K = 0):
    return np.exp( K*CosTheta )/(2*np.pi*sp.special.iv(0, K))


def von_Mises( K = 0):
    while True:
        CosTheta = 1-2*np.random.random()
        y = np.random.random()*von_Mises_PDF(1, K)
        if y <= von_Mises_PDF( CosTheta, K):
            return CosTheta

#K values
Kvalues = [0, 2, 10]
#Number of events
Nevents = 500

plt.rc('text', usetex=True)
plt.rc('font', size=16)


fig = plt.figure( figsize = (12,6) )
fig.subplots_adjust( hspace=0.0 )

for i in xrange(3):
    ax1 = fig.add_subplot(2,3,1+i,projection='hammer')
    ax2 = fig.add_subplot(2,1,2)#mollweide


    #Generating random sample and rotating
    Thetas = []
    Phis = []
    for i_e in xrange(Nevents):
        Phi = np.pi-2*np.pi*np.random.random()
        Phis.append( Phi )
        CosTheta = von_Mises(K = Kvalues[i])
        Thetas.append(np.pi/2-np.arccos(CosTheta))
        
    Thetas = np.array(Thetas)
    Phis = np.array(Phis)
    ax1.plot( Phis, Thetas, 'o', color='red', zorder = 10, markeredgewidth=0, ms=3 )

    ax1.set_xticks( np.linspace(-np.pi, np.pi, 6) )
    ax1.set_xticklabels( [''] )
    
    ax1.set_yticks( list(np.linspace(-75./90*np.pi/2, 75./90*np.pi/2, 5)))
    ax1.set_yticklabels( [165,135,90,45,15] )
    
    ax1.grid( linestyle = '-', linewidth=2 )

    costheta = np.linspace( -1,1,100 )
    ax2.plot( costheta, von_Mises_PDF( costheta, K = Kvalues[i]), '-', lw=3, label = '$k=%d$'%(Kvalues[i]) )
    ax2.legend( fancybox = True, loc = 'upper left' )
    #ax1.text( 0.7,0.8, '$k = %d$'%(Kvalues[i]), transform=ax2.transAxes, fontsize = 30 )
    ax1.set_title('$k=%d$'%(Kvalues[i]))
    ax2.set_yticks([])
    ax2.set_xticks(np.linspace(-1,1,5))
    #fig.gca().invert_xaxis()
    ax2.set_xticklabels(np.linspace(-1,1,5))
    if i != 2:
        ax2.set_xticklabels([])
    else:
        ax2.set_xlabel( '$\cos{\\theta}$', fontsize = 20 )
        ax2.set_ylabel( 'PDF', fontsize = 20 )
        #ax2.text( -0.1,2.1, 'PDF', transform=ax2.transAxes, fontsize = 14, rotation=90 )


#plt.show()
plt.savefig( "../vonMises.pdf" )
