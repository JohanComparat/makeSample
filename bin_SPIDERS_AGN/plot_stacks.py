import astropy.io.fits as fits
import time

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
import numpy as n

from cycler import cycler
# 1. Setting prop cycle on default rc parameter
p.rc('lines', linewidth=1.3)
p.rc('axes', prop_cycle=cycler('color', ["#638bd9", "#e586b6", "#dd7c25", "#598664", "#9631a9", "#cff159", "#534f55", "#ab1519", "#89dbde"]) )

import os
import sys

from scipy.stats import scoreatpercentile
from scipy.stats import norm
from scipy.interpolate import interp1d
m_bins = n.arange(-4,4,0.1)

import glob
path_2_spec_dir = os.path.join(os.environ['HOME'], 'SDSS/stacks/X_AGN')
fig_dir = os.path.join(os.environ['GIT_MAKESAMPLE'], 'figures/agn/figures_VAC')
#path_2_fly_dir =  os.path.join(os.environ['HOME'], 'SDSS/stacks/SPIDERS_C_GAL', 'firefly/')
file_list = n.array( glob.glob( os.path.join( path_2_spec_dir ,'*.stack')))
file_list.sort()
baseNames = n.array([ os.path.basename(bn)[5:-6] for bn in file_list ])
baseNamesSplit = n.array([ bn.split('_') for bn in baseNames ])
uniq_type = n.unique(baseNamesSplit.T[0])

def plot_stacks(file_list_i, baseNames_i, path_2_figure='x.png'):
	print(path_2_figure)
	fig = p.figure(0, (12.2, 8.2), frameon=False )
	# panel top left: spectrum + models
	fig.add_subplot(1,1,1, 
				title=baseNames_i[0].split('_')[0], 
				xlabel='wavelength [Angstrom, rest frame]', 
				ylabel=r'Flux [$f_\lambda$]', # 10^{-17}$ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]',
				#ylim=((n.min(y_data)*0.9,1.1* scoreatpercentile(y_data,90))),
				#xlim=((XMIN, XMAX)) 
				)
	p.grid()
	
	for jj, (path_2_spec, baseName) in enumerate(zip(file_list_i, baseNames_i)):
		d = fits.open(path_2_spec)
		sel = (d[1].data['NspectraPerPixel'] > 0.5*n.max(d[1].data['NspectraPerPixel'])) & (d[1].data['medianStack']>0)
		x_data = d[1].data['wavelength'][sel]
		y_data = d[1].data['medianStack'][sel]
		y_err = d[1].data['jackknifStackErrors'][sel]
		N_spec = d[1].data['NspectraPerPixel'][sel]
		labels = baseName.split('_')
		p.plot(x_data, y_data/2**jj, lw=1, label=str(float(labels[2])/10.)+'<z<'+str(float(labels[4])/10.) + ', N='+str(int(n.median(N_spec))) )

	p.yscale('log')
	p.legend(loc=0, fontsize=8, frameon=False)
	p.tight_layout()
	p.savefig(path_2_figure)
	p.clf()

for UT in uniq_type:
	sel = ( baseNamesSplit.T[0] == UT )
	plot_stacks(file_list[sel], baseNames[sel], path_2_figure=os.path.join(fig_dir, UT+'.png') )
	
	