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
#Data folder
DataFolder = '/home/bustamsn/PhD/DataHaswell/'
#DataFolder = '/home/bustamsn/PhD/DataFreya/'
#DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/'
#DataFolder = '/home/bustamsn/PhD/CosmologicalBlackHoles/Data'
DataResults = '.'
Resolution = 512
Simulation = 'cosmobh00/'

nshells = 60
zcut =  1                     # vertical cut in Kpc


#Morphology estimator       # Circularity    Fitting
estimator = 'Circularity'


#========================================================================================
#   GALAXY PLOT
#========================================================================================
attrs = ['pos', 'vel', 'mass', 'age', 'pot', 'id']
s = gadget_readsnap( 15, snappath='%s/Sims%3d/%s/output/snapdir_015'%(DataFolder, Resolution, Simulation), snapbase='snap_', hdf5=True, loadonlytype=[4,5], loadonly=attrs )
sf = load_subfind( 15, dir='%s/Sims%3d/%s/output/'%(DataFolder, Resolution, Simulation), hdf5=True, loadonly=['smty', 'SubhaloCM', 'flty', 'fnsh', 'slty', 'fpos', 'frc2', 'svel', 'spos', 'ffsh','fmty', 'ssph'] )
#Calculating indexes
IDs = np.arange(len(sf.data['fmty']))
mask_bh = sf.data['fmty'][:,5]>0
IDs_sub = np.insert(np.cumsum(sf.data['fnsh'])[:-1], 0, 0)

IDs = IDs[mask_bh]


IDs_sub = IDs_sub[mask_bh]
StellarMass2 = sf.data['fmty'][mask_bh][:,4]
StellarMass = sf.data['smty'][IDs_sub][:,4]
GasMass = sf.data['smty'][IDs_sub][:,0]
DMMass = sf.data['smty'][IDs_sub][:,1]
Colorgr = sf.data['ssph'][IDs_sub][:,4] - sf.data['ssph'][IDs_sub][:,5]
Spos = sf.data['spos'][IDs_sub]


mask_ord = np.argsort(StellarMass2)[::-1]

IDs = IDs[mask_ord]
IDs_sub = IDs_sub[mask_ord]
StellarMass = StellarMass[mask_ord]
GasMass = GasMass[mask_ord]
DMMass = DMMass[mask_ord]
Colorgr = Colorgr[mask_ord]
StellarMass2 = StellarMass2[mask_ord]
Spos = Spos[mask_ord]


Halos_params = []
for IDHALO in []:#xrange(len(IDs)):
    plt.close('all')
    s = gadget_readsnap( 15, snappath='%s/Sims%3d/%s/output/snapdir_015'%(DataFolder, Resolution, Simulation), snapbase='snap_', hdf5=True, loadonlytype=[4,5], loadonly=attrs )
    s.calc_sf_indizes( sf )
    s.select_halo( sf, use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True, haloid=IDs[IDHALO], subhalo=0 )

    g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][IDs[IDHALO]] )
    g.prep_data()

    sdata = g.sgdata['sdata']
    mass = sdata['mass']
    pos = sdata['pos']
    pot = sdata['pot']
    eps2 = sdata['eps2']


    #Using Circularities
    if estimator == 'Circularity':
        mask_bulge = g.sgdata['sdata']['eps2'] < 0
        Mass_Bulge = 2*np.sum(mass[mask_bulge])
        Mass_Disk = np.sum(mass) - Mass_Bulge
        Halos_params.append( Mass_Disk/(Mass_Disk+Mass_Bulge) )
        
        print IDHALO, Mass_Disk/(Mass_Disk+Mass_Bulge)



    #Using a Fitting of surface density profile (Exponential + Sersic)
    if estimator == 'Fitting':
        galrad = 0.1*sf.data['frc2'][IDs[IDHALO]]

        Rcut = galrad
        rd = np.linspace(0.0, Rcut, nshells)
        mnow = np.zeros(len(rd))

        rad = pylab.sqrt( (pos[:,1:]**2).sum(axis=1) )
        z = pos[:,0]

        ii, = np.where( (abs(z) < zcut) )
        weight = mass*1e10/0.7

        bins = nshells
        sden, edges = np.histogram( rad[ii], bins=bins, range=(0., Rcut), weights=weight[ii] )
        sa = np.zeros(len(edges)-1)
        sa[:] = np.pi * (edges[1:]**2 - edges[:-1]**2)
        sden /= sa
        sden *= 1e-6
        x = np.zeros(len(edges)-1)
        x[:] = 0.5 * (edges[1:] + edges[:-1])
        r = x

        bounds = ([0., 0., 0., 0., 0.1], [1e4, 30., 1e4, 10, 5.])
        sigma = 0.1*sden
        guess = [1e2,10,1e2,1,1]
        #guess = [1,10,1,1,1]
        #(params, pcov) = curve_fit( Total_Profile, r, sden, guess, sigma=sigma, bounds=bounds )
        #(params, pcov) = curve_fit( Total_Profile, r, sden, guess, bounds=bounds )
        #(params, pcov) = curve_fit( Total_Profile, r, sden, guess, sigma=sigma )
        #(params, pcov) = curve_fit( Total_Profile, r, sden, guess )
        
        try:
            (params, pcov) = curve_fit( Total_Profile, r, sden, guess, sigma=sigma, bounds=bounds )
        except:
            try:
                (params, pcov) = curve_fit( Total_Profile, r, sden, guess )
            except:
                params = guess
                
        
        Halos_params.append( params )

        print params

        rarray = np.linspace( 0, Rcut, 500 )

        plt.semilogy( r, sden, 'o', lw=2, color='black' )
        #plt.semilogy( rarray, total_profile(rarray, params[0], params[1], params[2], params[3], params[4]), color='red', lw=2 )
        #plt.semilogy( rarray, exp_prof(rarray, params[0], params[1]), color='blue', lw=2, label='Disk' )
        #plt.semilogy( rarray, sersic_prof1(rarray, params[2], params[3], params[4]), color='green', lw=2, label='Bulge' )

        plt.semilogy( rarray, Total_Profile(rarray, params[0], params[1], params[2], params[3], params[4]), color='red', lw=2 )
        plt.semilogy( rarray, Exponential_Disk(rarray, params[0], params[1]), color='blue', lw=2, label='Disk' )
        plt.semilogy( rarray, Sersic_Profile(rarray, params[2], params[3], params[4]), color='green', lw=2, label='Bulge' )

        plt.ylim( (1e-2, 1e5) )
        plt.legend()
        plt.savefig( './figs/morphologies/Halo_%d.png'%IDHALO )
    

#Halos_params = np.array(Halos_params)

#np.savetxt( './data/morphologies_cosmobh04.txt', Halos_params )



#DT = abs(np.loadtxt( './data/morphologies_cosmobh04.txt' ))


#plt.figure( figsize = (14,6) )

#plt.subplot( 1,2,1 )
#plt.scatter( np.log10(StellarMass) + 10 , Colorgr, alpha = 0.5, marker = 'o', s=20, c=DT, vmin=0, vmax=1 )
#plt.xlim( [8.5,12] )
#plt.ylim( [0.0,1.0] )
#cb1 = plt.colorbar()
#plt.xlabel( '$\log M_{\star} [10^{10}\ M_{\odot}]$', fontsize = 18 )
#plt.ylabel( 'g - r', fontsize = 18 )
#cb1.set_label( 'D / T', fontsize = 18 )


#plt.subplot( 1,2,2 )
#plt.scatter( np.log10(StellarMass2) + 10 , DT, c=Colorgr, vmin = 0.0, vmax = 1.0, alpha = 0.5, marker = 'o', s=20)
#cb2 = plt.colorbar()
#plt.ylim( [0.0,1.0] )
#plt.xlim( [8.5,12] )
#plt.xlabel( '$\log M_{\star} [10^{10}\ M_{\odot}]$', fontsize = 18 )
#plt.ylabel( 'D / T', fontsize = 18 )
#cb2.set_label( 'g - r', fontsize = 18 )

##plt.savegig('./fig.pdf')

#plt.show()

#3D plot

#fig = plt.figure()
#ax = fig.gca(projection='3d')

#mask = g.sgdata['sdata']['eps2'] > 0.7

#ax.plot( pos[mask][:,1], pos[mask][:,2], pos[mask][:,0], 'o', ms=2, color='red' )
#ax.plot( pos[mask==False][:,1], pos[mask==False][:,2], pos[mask==False][:,0], '.', ms=1 )

#ax.set_aspect(1)
#plt.show()
