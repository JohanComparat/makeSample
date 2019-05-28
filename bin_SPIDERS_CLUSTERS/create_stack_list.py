
import sys
import numpy as n
from scipy.interpolate import interp1d
import os
import astropy.io.fits as fits
import pickle

version = 'v1'

zmins = n.hstack((0., 0.005, n.arange(0.1,3.0,0.1) ))
zmaxs = n.hstack((0.005, n.arange(0.,3.0,0.1)+0.5))

# on laptop
#path_2_cats = os.path.join(os.environ['HOME'], 'data', 'spiders', 'agn')
# on servers
path_2_cats = os.path.join(os.environ['HOME'], 'hegcl/SPIDERS')
path_2_stack_lists = os.path.join(os.environ['HOME'], 'SDSS/stacks/SPIDERS_C_GAL')

# 4 catalogs 
path_2_CAT = os.path.join(path_2_cats, 'validatedclusters_catalogue_2018-04-27_version_round123-v1_Xmass123-v1-flat.fits')

cat = fits.open(path_2_CAT)[1].data

import astropy.cosmology as co
aa=co.Planck15
import astropy.units as u
from astropy.coordinates import angles
#import AngularSeparation
from astropy import coordinates as coord

center = coord.ICRS(ra=cat['RA_OPT']*u.degree, dec=cat['DEC_OPT']*u.degree)
position = coord.ICRS(cat['RA_GAL']*u.degree, cat['DEC_GAL']*u.degree)
sep_r200c =  (center.separation(position)/(cat['R200C_DEG']*u.degree)).value 

z_clus = cat['CLUZSPEC']

ok = (sep_r200c>=0) & (cat['ISMEMBER']==1) & (0.1<z_clus) & (z_clus<0.3)

delta_z = abs((cat['IDLSPEC1D_Z_NOQSO'] - z_clus)*3e5/cat['CLUVDISPBEST'] )


log_sep_r200c = n.log10(sep_r200c[ok])
bins = n.arange(n.min(log_sep_r200c)-1, n.max(log_sep_r200c)+2, 0.00001)

out_N = n.cumsum(n.histogram(log_sep_r200c, bins=bins)[0])
from scipy.interpolate import interp1d
itp = interp1d( out_N, bins[1:] )
	
bins_x = n.hstack((itp( n.arange(0, len(log_sep_r200c), 400 ) ), log_sep_r200c.max()+0.1))

for bin_low, bin_high in zip (bins_x[:-1], bins_x[1:]):
	name = "clusterCGAL_"+str(bin_low)+"_r_"+str(bin_high)+'.ascii'
	print(name)
	selection = (n.log10(sep_r200c)>=bin_low) & (n.log10(sep_r200c)<bin_high) & (sep_r200c>=0) & (cat['ISMEMBER']==1)
	DATA = n.transpose([cat['PLATE'], cat['MJD'], cat['FIBERID'], cat['IDLSPEC1D_Z_NOQSO']])[selection]
	if len(cat['PLATE'][selection])>10:
		n.savetxt(os.path.join(path_2_stack_lists, name), DATA)


for bin_low, bin_high in zip (bins_x[:-1], bins_x[1:]):
	name = "clusterCGAL_"+str(bin_low)+"_r_"+str(bin_high)+'_lowCZ.ascii'
	print(name)
	selection = (n.log10(sep_r200c)>=bin_low) & (n.log10(sep_r200c)<bin_high) & (sep_r200c>=0) & (cat['ISMEMBER']==1)&(delta_z<1./3.)
	DATA = n.transpose([cat['PLATE'], cat['MJD'], cat['FIBERID'], cat['IDLSPEC1D_Z_NOQSO']])[selection]
	if len(cat['PLATE'][selection])>10:
		n.savetxt(os.path.join(path_2_stack_lists, name), DATA)

for bin_low, bin_high in zip (bins_x[:-1], bins_x[1:]):
	name = "clusterCGAL_"+str(bin_low)+"_r_"+str(bin_high)+'_highCZ.ascii'
	print(name)
	selection = (n.log10(sep_r200c)>=bin_low) & (n.log10(sep_r200c)<bin_high) & (sep_r200c>=0) & (cat['ISMEMBER']==1)&(delta_z<1./3.)
	DATA = n.transpose([cat['PLATE'], cat['MJD'], cat['FIBERID'], cat['IDLSPEC1D_Z_NOQSO']])[selection]
	if len(cat['PLATE'][selection])>10:
		n.savetxt(os.path.join(path_2_stack_lists, name), DATA)

