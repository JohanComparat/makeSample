import astropy.io.fits as fits
import os, sys, glob
from os.path import join
#import pymangle as mangle
import numpy as np
import matplotlib.pyplot as p
from scipy.interpolate import interp1d

import astropy.units as u
import astropy.cosmology as cc
cosmo = cc.Planck13
from astropy.coordinates import SkyCoord

from sklearn.neighbors import BallTree, DistanceMetric
from astropy.table import Table,unique
from math import radians, cos, sin, asin, sqrt, pi

in_dir  = '/data36s/comparat/CODEX_clustering/catalogs/'
out_dir = '/data36s/comparat/CODEX_clustering/angular_clustering/'

deg_to_rad = np.pi/180.
arcsec = 1. / 3600.
rs = 10**np.arange(-1,1.6,0.1) *arcsec
#rs = 10**np.arange(-1,2.6,0.1) *arcsec

# DATA 
path_2_data_2rxs   = join( in_dir, 'cat_spiders_masked_X.fits' )
data_file = path_2_data_2rxs

ra_name_data = 'RA'
dec_name_data = 'Dec'

hduD     = fits.open(data_file)
ra_data    = hduD[1].data[ra_name_data]
dec_data    = hduD[1].data[dec_name_data]
z_data  = hduD[1].data['z']

coords = SkyCoord(ra_data, dec_data, unit='deg', frame='icrs')
bb_data = coords.galactic.b.value
ll_data = coords.galactic.l.value
bb_ecl_data = coords.barycentrictrueecliptic.lat

selection_data = (z_data>0)

N_data = len(ra_data[selection_data])

# GAIA CATALOGS
gaia_dir = '/data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/GAIA/DR2/'
gaia_table_list = np.array(glob.glob(os.path.join(gaia_dir, 'table_*.fits')))
gaia_table_list.sort()

for gaia_file in gaia_table_list[1:][::-1][:5]:
	print(gaia_file)
	hdu_G     = fits.open(gaia_file)
	ra_gaia, dec_gaia = hdu_G[1].data['ra'], hdu_G[1].data['dec']
	coords = SkyCoord(ra_gaia, dec_gaia, unit='deg', frame='icrs')
	bb_gaia = coords.galactic.b.value
	ll_gaia = coords.galactic.l.value
	bb_ecl_gaia = coords.barycentrictrueecliptic.lat

	x_gal_gaia = (abs(bb_gaia)>20)&(dec_gaia<80)&(dec_gaia>-80)&(bb_ecl_gaia.value>-80)
	selection_gaia = (x_gal_gaia)
	N_gaia = len(ra_gaia[selection_gaia])


	agn_coordinates = deg_to_rad * np.transpose([dec_data[selection_data], ra_data[selection_data]])
	gaia_coordinates = deg_to_rad * np.transpose([dec_gaia[selection_gaia], ra_gaia[selection_gaia]])

	Tree_obj_Gaia = BallTree(gaia_coordinates, metric='haversine') 
	Tree_obj_AGN = BallTree(agn_coordinates , metric='haversine') 

	test_c = np.array([ Tree_obj_Gaia.query_radius(agn_coordinates, r = rr, count_only=True) for rr in rs ])

	N_pairs_total = test_c.sum(axis=1)
	Delta_N_pairs = N_pairs_total[1:]-N_pairs_total[:-1]
	area = 4.*np.pi*(rs[1:]**2 - rs[:-1]**2)

	pair_density = Delta_N_pairs/(area*N_data*N_gaia)

	out_data = os.path.join(out_dir , 'cat_spiders_masked_X_GAIA_'+os.path.basename(gaia_file)+'.data')
	np.savetxt(out_data, np.transpose([rs[1:], pair_density]), header='radius_arcsec density' )



# EXCEPTION for the first file, that has a broad magnitude range 1-5 that I re-cut by hand into 2 and extend the radius.

gaia_file = gaia_table_list[0]
print(gaia_file)
hdu_G     = fits.open(gaia_file)
ra_gaia, dec_gaia = hdu_G[1].data['ra'], hdu_G[1].data['dec']
coords = SkyCoord(ra_gaia, dec_gaia, unit='deg', frame='icrs')
bb_gaia = coords.galactic.b.value
ll_gaia = coords.galactic.l.value
bb_ecl_gaia = coords.barycentrictrueecliptic.lat
g_mag = hdu_G[1].data['phot_g_mean_mag']

x_gal_gaia = (abs(bb_gaia)>20)&(dec_gaia<80)&(dec_gaia>-80)&(bb_ecl_gaia.value>-80)


mag_sel = (g_mag>4)
selection_gaia = (x_gal_gaia)&(mag_sel)
N_gaia = len(ra_gaia[selection_gaia])

# Tree to select pairs

# COUNT UNIQUE 

agn_coordinates = deg_to_rad * np.transpose([dec_data[selection_data], ra_data[selection_data]])
gaia_coordinates = deg_to_rad * np.transpose([dec_gaia[selection_gaia], ra_gaia[selection_gaia]])

Tree_obj_Gaia = BallTree(gaia_coordinates, metric='haversine') 
Tree_obj_AGN = BallTree(agn_coordinates , metric='haversine') 

test_c = np.array([ Tree_obj_Gaia.query_radius(agn_coordinates, r = rr, count_only=True) for rr in rs ])

N_pairs_total = test_c.sum(axis=1)
Delta_N_pairs = N_pairs_total[1:]-N_pairs_total[:-1]
area = 4.*np.pi*(rs[1:]**2 - rs[:-1]**2)

pair_density = Delta_N_pairs/(area*N_data*N_gaia)

out_data = os.path.join(out_dir , 'cat_spiders_masked_X_GAIA_table_4_g_5.fits.data')
np.savetxt(out_data, np.transpose([rs[1:], pair_density]), header='radius_arcsec density' )

mag_sel = (g_mag<=4)
selection_gaia = (x_gal_gaia)&(mag_sel)
N_gaia = len(ra_gaia[selection_gaia])

# Tree to select pairs

# COUNT UNIQUE 

agn_coordinates = deg_to_rad * np.transpose([dec_data[selection_data], ra_data[selection_data]])
gaia_coordinates = deg_to_rad * np.transpose([dec_gaia[selection_gaia], ra_gaia[selection_gaia]])

Tree_obj_Gaia = BallTree(gaia_coordinates, metric='haversine') 
Tree_obj_AGN = BallTree(agn_coordinates , metric='haversine') 

test_c = np.array([ Tree_obj_Gaia.query_radius(agn_coordinates, r = rr, count_only=True) for rr in rs ])

N_pairs_total = test_c.sum(axis=1)
Delta_N_pairs = N_pairs_total[1:]-N_pairs_total[:-1]
area = 4.*np.pi*(rs[1:]**2 - rs[:-1]**2)

pair_density = Delta_N_pairs/(area*N_data*N_gaia)

out_data = os.path.join(out_dir , 'cat_spiders_masked_X_GAIA_table_1_g_4.fits.data')
np.savetxt(out_data, np.transpose([rs[1:], pair_density]), header='radius_arcsec density' )

