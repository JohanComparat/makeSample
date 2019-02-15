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
from astropy.table import Table,unique,Column
from math import radians, cos, sin, asin, sqrt, pi

in_dir  = '/data36s/comparat/AGN_clustering/catalogs/'
out_dir = '/data36s/comparat/AGN_clustering/angular_clustering/'

deg_to_rad = np.pi/180.
arcsec = 1. / 3600.
rs = np.arange(0,20,.1)*arcsec

# SDSS WISE CATALOGS
# path_2_file = '/data44s/eroAGN_WG_DATA/DATA/masks/SDSS_WISE_imageprop_nside512.fits'
# hdu_S     = fits.open(path_2_file)

# DATA 
path_2_data_2rxs   = join( in_dir, '2RXS_AllWISE_catalog_paper_2017May26_v5_11_0_sdss_26_VERON_MASKED.fits'   )

out_file = path_2_data_2rxs[:-5]+'_GAIA_star_mask.fits'

data_file = path_2_data_2rxs

ra_name_data = 'ALLW_RA'
dec_name_data = 'ALLW_DEC'

hduD = Table.read(data_file)

ra_data    = hduD[ra_name_data]
dec_data    = hduD[dec_name_data]

agn_coordinates = deg_to_rad * np.transpose([dec_data, ra_data])
Tree_obj_AGN = BallTree(agn_coordinates , metric='haversine') 

# GAIA CATALOGS
gaia_dir = '/data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/GAIA/DR2/'
gaia_table_list = np.array(glob.glob(os.path.join(gaia_dir, 'table_*.fits')))
gaia_table_list.sort()

valid_ids = np.array([0, 1, 2,3,4,8,9,10,11,12])

gaia_file_2_radius = {
	'table_-10_g_5.fits' :  10. ,
	'table_5_g_6.fits'   :  6 ,
	'table_6_g_7.fits'   :  4 ,
	'table_7_g_8.fits'   :  3 ,
	'table_8_g_9.fits'   :  2 ,
	'table_9_g_10.fits'  :  2 ,
	'table_10_g_11.fits' :  2 ,
	'table_11_g_12.fits' :  2 ,
	'table_12_g_13.fits' :  2 ,
	'table_13_g_14.fits' :  2 ,
	'table_14_g_15.fits' :  1 ,
	'table_15_g_16.fits' :  1 ,
	'table_16_g_17.fits' :  1 
	}

star_nearby = []

for gaia_file in gaia_table_list[valid_ids]:
	print(gaia_file)
	r_mask = gaia_file_2_radius[os.path.basename(gaia_file)] / 3600.
	hdu_G     = fits.open(gaia_file)
	ra_gaia, dec_gaia = hdu_G[1].data['ra'], hdu_G[1].data['dec']
	gaia_coordinates = deg_to_rad * np.transpose([dec_gaia, ra_gaia])
	print('measures distances')
	Tree_obj_Gaia = BallTree(gaia_coordinates, metric='haversine') 
	test_c = Tree_obj_Gaia.query_radius(agn_coordinates, r = r_mask, count_only = True) 
	to_be_masked = (test_c>0)
	star_nearby.append( to_be_masked )
	print('N to mask:', len(to_be_masked.nonzero()[0]))

mask_array = np.transpose(star_nearby).any(axis=1)
c_mask = Column(mask_array, name='nearby_gaia_star')

hduD.add_column(c_mask)

print( out_file )
if os.path.isfile(out_file):
	os.system("rm "+out_file)

hduD.write(out_file, format='fits')
