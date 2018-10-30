from pylab import *
from loadmodules import *
from const import *
from util import *
from parse_particledata import parse_particledata
from scipy.optimize import curve_fit
from lmfit import Model
from scipy.special import gamma
import sys
p = plot_helper.plot_helper()

colors = ['g','b','r','k','c','y','m','purple']
lws = [0.7, 1., 1.3]

band    = {'Umag':0, 'Bmag':1, 'Vmag':2, 'Kmag':3, 'gmag':4, 'rmag':5, 'imag':6, 'zmag':7}
Msunabs = [5.61, 5.48, 4.83, 3.30, 5.12, 4.68, 4.57, 4.60]

def plot_stellar_surfden( runs, dirs, outpath, outputlistfile, zlist, suffix, nrows, ncols, iband=None, writepfile=False, restest=False, readguessparams=True, allowonecompfit=False):

	if len(zlist) == 1:
		if restest:
			panels = len(runs) / len(restest)
		else:
			panels = nrows*ncols

		figure = multipanel_layout.multipanel_layout( nrows=nrows, ncols=ncols, npanels=panels )
		figure.set_fontsize(3)
		figure.set_figure_layout()
		figure.set_axis_locators( xminloc=5.0, xmajloc=5.0, yminloc=0.1, ymajloc=0.5 )
		figure.set_fontsize(10)
		figure.set_axis_labels( xlabel="$\\rm{R [kpc]}$", ylabel="$\\rm{\Sigma [M_{\odot} pc^{-2}]}$" )
		figure.set_axis_limits_and_aspect( xlim=[0.0, 19.0], ylim=[1e0, 2e4], logyaxis=True )

	clist = ['r', 'b', 'c', 'g', 'k']

	apath = outpath #+ '/plotall/'

	if writepfile and len(zlist) == 1:
		filename = apath + '/fit_table.txt'
		f2 = open(filename, 'w')
		header = "%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n" % ("Run", "Time", "virialmass", "virialradius", "Stellarmass", "discmass", "Rd", "bulgemass", "reff", "n", "D/T")
		f2.write(header)

	#filename = outpath + '/r_optical.txt'
	#data = np.loadtxt(filename, comments='#')
	#hname = data[:,0]
	#roptmin = data[:,1]
	#roptmax = data[:,2]
	hname = np.array(runs)
	roptmin = [30.] * len(runs)
	#roptmin = np.array([30., 33., 30., 30., 32., 30.])

	if readguessparams:
		filename = outpath + '/discbulgefitparam_initialguess.txt'
		gudat = np.loadtxt(filename, comments='#', dtype={'names': ('halo','gp1','gp2','gp3','gp4','gp5','sdlim'), 'formats': ('S3','f4','f4','f4','f4','f4','f4')})
		print "gudat=",gudat
		hname = gudat['halo']
		print "hname=",hname
		gparams = np.stack((gudat['gp1'],gudat['gp2'],gudat['gp3'],gudat['gp4'],gudat['gp5']), axis=1)
		gsdlim = gudat['sdlim']
		print "gsdlim=",gsdlim
		print "gparams=",gparams
		

	p = plot_helper.plot_helper()

        for d in range( len( runs ) ): 
		if len(zlist) > 1:
			panels = nrows*ncols

			figure = multipanel_layout.multipanel_layout( nrows=nrows, ncols=ncols, npanels=panels )
			figure.set_figure_layout()
			figure.set_axis_locators( xminloc=5.0, xmajloc=5.0 )
			figure.set_axis_labels( xlabel="$\\rm{R [kpc]}$", ylabel="$\\rm{\Sigma [M_{\odot} pc^{-2}]}$" )
			figure.set_axis_limits_and_aspect( xlim=[0.0, 19.0], ylim=[1e0, 2e4], logyaxis=True )
			figure.set_fontsize()

                dd = dirs[d] + runs[d]

		if not restest:
			ax = figure.axes[d]
		elif len(zlist) == 1:
			ax = figure.axes[int(d/len(restest))]

		wpath = outpath+runs[d]+'/radprof/'
		if not os.path.exists(wpath):
			os.makedirs(wpath)

		if writepfile and len(zlist) > 1:
			filename = wpath + 'fit_table_%s.txt' % runs[d]
			f2 = open(filename,'w')
			header = "%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n" % ("Run", "Time", "virial mass", "virial radius", "Stellar mass", "disc mass", "R_d", "bulge mass", "r_eff", "n", "D/T", "R_dg", "gagdisc mass")
			f2.write(header)

		nshells = 60 # 35 up to galrad is OK

		print "outputlistfile=",outputlistfile
		print "zlist=",zlist
		snaps = select_snapshot_number.select_snapshot_number( outputlistfile[d], zlist )
		if isinstance(snaps, int):
                        snaps = [snaps]

		#parameters = np.linspace(0.,8.,len(snaps))
		#cmap = mpl.cm.jet
		#s_m = figure.get_scalar_mappable_for_colorbar(parameters, cmap)
		time = np.zeros(len(snaps))

		for isnap in range(len(snaps)):
			if len(zlist) > 1:
				ax = figure.axes[isnap]

			snap = snaps[isnap]
			print "Doing dir %s snap %d." % (dd, snap)

			attrs = ['pos', 'vel', 'mass', 'age']
			if iband:
				attrs.append('gsph')
			s = gadget_readsnap( snap, snappath=dd+'/output/', hdf5=True, loadonlytype=[4], loadonly=attrs )
			sf = load_subfind( snap, dir=dd+'/output/', hdf5=True, loadonly=['fpos', 'frc2', 'svel', 'flty', 'fnsh', 'slty', 'fmc2'] )
			s.calc_sf_indizes( sf )
			s.select_halo( sf, [], 3., use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True )
			
			#attrs = ['pos', 'vel', 'mass', 'age', 'id']
			g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][0] )
			g.prep_data()

			sdata = g.sgdata['sdata']
			mass = sdata['mass']
			pos = sdata['pos']

			galrad = 0.1 * sf.data['frc2'][0]
			time[isnap] = s.cosmology_get_lookback_time_from_a( s.time.astype('float64'), is_flat=True )
			if time[isnap] < 0.0001:
				time[isnap] = 0.0
			t0 = str(time[isnap])
			t0l = 'ti='+t0[:4]+' Gyr'

			zcut =  0.001                      # vertical cut in Mpc
			Rcut =  0.040 

                        rd = np.linspace(0.0, Rcut, nshells)
                        mnow = np.zeros(len(rd))

			rad = pylab.sqrt( (pos[:,1:]**2).sum(axis=1) )
			z = pos[:,0]

			ii, = np.where( (abs(z) < zcut) )

			if iband:
				luminosity = 10**( -0.4 * (sdata[iband] - Msunabs[band[iband]]) )
				weight = luminosity
			else:
				weight = mass

			#bins = np.logspace(0., np.log10(Rcut*1e3), nshells) * 1e-3
			#print "bins=",bins
			bins = nshells
			sden, edges = np.histogram( rad[ii], bins=bins, range=(0., Rcut), weights=weight[ii] )
			sa = np.zeros(len(edges)-1)
			sa[:] = np.pi * (edges[1:]**2 - edges[:-1]**2)
			sden /= sa

			x = np.zeros(len(edges)-1)
			x[:] = 0.5 * (edges[1:] + edges[:-1])

			sden *= 1e-6
			r = x * 1e3

                        print "Surface density at R=8 kpc:",sden[find_nearest(r, [8.])]

			ropt = roptmin[hname==runs[d]]
			print "galrad,ropt=",galrad,ropt
			indy = find_nearest( x, [ropt*1e-3] ).astype('int64')	
			if runs[d] == 'halo_11': # special case: companion
				indy = find_nearest( x, [0.07*sf.data['frc2'][0]] ).astype('int64')
			else:
                                indy = find_nearest( x, [galrad] ).astype('int64')
				
			if not readguessparams:
				sdlim = 1.
			else:
				ind = np.where( hname == runs[d].split('_')[1] )
				sdlim = gsdlim[ind].ravel()[0]
				print "sdlim=",sdlim
			indy = find_nearest( sden*1e4, [sdlim] ).astype('int64')

			rfit = x[indy]*1e3
			sdfit = sden[:indy]
			r = r[:indy][sdfit>0.]
			sdfit = sdfit[sdfit>0.]

			print "fitting carried out in radial range %f - %f" % (x[0]*1e3,x[indy]*1e3)

			if not restest:
				flag = 0
				try:
					if readguessparams:
						ind = np.where( hname == runs[d].split('_')[1] )
						guess = gparams[ind,:].ravel()
					# guess params: (exp amp, Rd, sersic amp, sersic b, 1/n)
					else:
						guess = (1., 4., 1., 1., 1.)#, 0.1, 10.)

					bounds = ([0., 0., 0., 0., 0.1], [10., 30., 10., 10., 5.])
					sigma = 0.1*sdfit
					(popt, pcov) = curve_fit( p.total_profile, r, sdfit, guess, sigma=sigma, bounds=bounds )

					# compute total masses from the fit                                                     
					disk_mass = 2.0 * np.pi * popt[0] * popt[1] * popt[1]
					bulge_mass = np.pi * popt[2] * popt[3] * popt[3] * gamma(2.0 /popt[4] + 1)
					disk_to_total = disk_mass / (bulge_mass + disk_mass)
					# Set to nan any crazy values (can happen at earlier times)              
					plist = np.array([disk_mass, bulge_mass, popt[0], popt[1], popt[2], popt[3], popt[4], disk_to_total])
					aa, = np.where( (plist > 30.) | (plist < 0.) )
					plist[aa] = nan
					disk_mass, bulge_mass, popt[0], popt[1], popt[2], popt[3], popt[4], disk_to_total = plist[0], plist[1], plist[2], plist[3], plist[4], plist[5], plist[6], plist[7]
					param = "%12s%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.2f\n" % (("Au%s"%runs[d].split('o')[1]), time[isnap], sf.data['fmc2'][0], sf.data['frc2'][0]*1e3, mass.sum(), disk_mass, popt[1], bulge_mass, popt[3] * p.sersic_b_param(1.0 / popt[4])**(1.0 / popt[4]), 1.0 / popt[4], disk_to_total)#, poptg[1], gasdisk_mass)
					if writepfile:
						f2.write(param)
					perr = np.sqrt(np.diag(pcov)).mean()
					if allowonecompfit:

						try:
							bguess = guess[2:]
							bounds = ([0., 0., 0.5], [10., 10., 4.])
							(popt2, pcov2) = curve_fit( p.sersic_prof1, r, sdfit, bguess, sigma=sigma, bounds=bounds )
							perr2 = np.sqrt(np.diag(pcov2)).mean()
							if perr2 < perr:
								print "**bulge only profile preferred"
								flag = 1
								popt = popt2
								pcov = pcov2
								perr = perr2
								disk_mass = nan
								bulge_mass = np.pi * popt[0] * popt[1] * popt[1] * gamma(2.0 /popt[2] + 1)
								disk_to_total = nan
						except:
							pass

					# disc only
						try:
							dguess = guess[:2]
							bounds = ([0., 0.], [10., 30.])
							(popt2, pcov2) = curve_fit( p.exp_prof, r, sdfit, dguess, sigma=sigma, bounds=bounds )
							perr2 = np.sqrt(np.diag(pcov2)).mean()
							if perr2 < perr:
								print "**disc only profile preferred"
								flag = 2
								popt = popt2
								pcov = pcov2
								perr = perr2
								disk_mass = 2.0 * np.pi * popt[0] * popt[1] * popt[1]
								bulge_mass = nan
								disk_to_total = nan
						except:
							pass
	
					print "popt=",popt
					print "pcov=",pcov
					print "perr=",perr
				except:
					popt = np.zeros(5)
					disk_mass = nan
					bulge_mass = nan
					disk_to_total = nan
					print "fitting failed"

				indy = find_nearest( x, rd ).astype('int64')
				nrad = len(rd)
				rm = x[indy]
				mnow[:] = sden[indy]
				mprev = mnow

				ax.semilogy(r,1e10*sdfit * 1e-6,'o',markersize=1.5,color='k', linewidth=0.)
				if flag == 1:
					ax.semilogy(r,1e10 * p.sersic_prof1(r, popt[0], popt[1], popt[2]) * 1e-6,'r--', label=r'n=%.2f' % (1./popt[2]))
				if flag == 0 or flag == 2:
					ax.semilogy(r,1e10 * p.exp_prof(r, popt[0], popt[1]) * 1e-6,'b-',label=r'$R_d = %.2f$' % (popt[1]))
				if flag == 0:
					ax.semilogy(r,1e10 * p.sersic_prof1(r, popt[2], popt[3], popt[4]) * 1e-6,'r--', label=r'n=%.2f' % (1./popt[4]))
					ax.semilogy(r,1e10*p.total_profile(r, popt[0], popt[1], popt[2], popt[3], popt[4]) * 1e-6, 'k-')

				ax.axvline(rfit, color='gray', linestyle='--')
				ax.xaxis.set_ticks([0,5,10,15,20])
				ax.legend( loc='upper right', frameon=False, prop={'size':6} )
				if len(zlist) > 1:
					figure.set_panel_title( panel=isnap, title="$\\rm %s\, %.2f$" % ("Time=",time[isnap]), position='bottom right') 
				else:
					figure.set_panel_title(panel=d, title="$\\rm{%s\,{%s}}$" % ("Au", runs[d].split('_')[1]), position='bottom left', color='k', fontsize=8 )
			else:
				for i, level in enumerate(restest):
					lst = 'level%01d' % level
					if lst in dirs[d]:
						label = " $\\rm{Au %s lvl %01d}$" % (runs[d].split('_')[1], level )
						lw = lws[i]
				ax.semilogy(r,1e10*sdfit * 1e-6, '-', lw=lw, color=colors[d%len(colors)], label=label)

                ax.legend( loc='upper left', frameon=False, prop={'size':7} )

		if len(zlist) > 1:
			print "filename=%s/starden_profevol_%s.%s"%(wpath, runs[d], suffix)
			figure.fig.savefig( '%s/starden_profevol_%s.%s' % (wpath, runs[d], suffix) )
			if writepfile:
				f2.close()

			#rtime = [ round(elem, 1) for elem in time ]

	if len(zlist) == 1:
		figure.reset_axis_limits()	
		print "Saving figures"
		if not restest:
			figure.fig.savefig( '%s/starden_prof%03d.%s' % (apath, snap, suffix) )
		else:
			figure.fig.savefig( '%s/starden_prof%03d_restest.%s' % (apath, snap, suffix) )
		if writepfile:
			f2.close()

def find_nearest(array,value):
        if len(value) == 1:
                idx = (np.abs(array-value)).argmin()
        else:
		idx = np.zeros(len(value))
                for i in range(len(value)):
                        idx[i] = (np.abs(array-value[i])).argmin()

        return idx
