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

in_dir  = '/data36s/comparat/AGN_clustering/catalogs/'
out_dir = '/data36s/comparat/AGN_clustering/angular_clustering/'


# SDSS WISE CATALOGS
# path_2_file = '/data44s/eroAGN_WG_DATA/DATA/masks/SDSS_WISE_imageprop_nside512.fits'
# hdu_S     = fits.open(path_2_file)

# DATA 
path_2_data_2rxs   = join( in_dir, '2RXS_AllWISE_catalog_paper_2017May26_v5_11_0_sdss_26_VERON_MASKED.fits'   )
data_file = path_2_data_2rxs

ra_name_data = 'ALLW_RA'
dec_name_data = 'ALLW_DEC'

hduD     = fits.open(data_file)
ra_data    = hduD[1].data[ra_name_data]
dec_data    = hduD[1].data[dec_name_data]
z_data = np.zeros_like(ra_data)
ratelim_data    = hduD[1].data['RATELIM']

coords = SkyCoord(ra_data, dec_data, unit='deg', frame='icrs')
bb_data = coords.galactic.b.value
ll_data = coords.galactic.l.value
bb_ecl_data = coords.barycentrictrueecliptic.lat

stars_data = (hduD[1].data['p_any']>0.5)&(hduD[1].data['2RXS_ExiML']>10)
rt_sel_data = (ratelim_data>0.015)
x_gal_data = (abs(bb_data)>20)&(dec_data<80)&(dec_data>-80)&(bb_ecl_data.value>-80)
selection_data = (x_gal_data)&(stars_data==False)&(rt_sel_data)

N_data = len(ra_data[selection_data])

# GAIA CATALOGS
gaia_dir = '/data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/GAIA/DR2/'
gaia_table_list = np.array(glob.glob(os.path.join(gaia_dir, 'table_*.fits')))
gaia_table_list.sort()

for gaia_file in gaia_table_list[1:]:
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

	# Tree to select pairs

	# COUNT UNIQUE 

	from sklearn.neighbors import BallTree, DistanceMetric
	from astropy.table import Table,unique
	from math import radians, cos, sin, asin, sqrt, pi

	deg_to_rad = np.pi/180.

	arcsec = 1. / 3600.
	rs = np.arange(100)*arcsec

	agn_coordinates = deg_to_rad * np.transpose([dec_data[selection_data], ra_data[selection_data]])
	gaia_coordinates = deg_to_rad * np.transpose([dec_gaia[selection_gaia], ra_gaia[selection_gaia]])

	Tree_obj_Gaia = BallTree(gaia_coordinates, metric='haversine') 
	Tree_obj_AGN = BallTree(agn_coordinates , metric='haversine') 

	test_c = np.array([ Tree_obj_Gaia.query_radius(agn_coordinates, r = rr, count_only=True) for rr in rs ])

	N_pairs_total = test_c.sum(axis=1)
	Delta_N_pairs = N_pairs_total[1:]-N_pairs_total[:-1]
	area = 4.*np.pi*(rs[1:]**2 - rs[:-1]**2)

	pair_density = Delta_N_pairs/(area*N_data*N_gaia)

	out_data = os.path.join(out_dir , '2RXS_AllWISE_catalog_paper_2017May26_X_GAIA_'+os.path.basename(gaia_file)+'.data')
	np.savetxt(out_data, np.transpose([rs[1:], pair_density]), header='radius_arcsec density' )

