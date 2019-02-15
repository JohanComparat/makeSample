import astropy.io.fits as fits
import os, sys
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

# catalogs :
path_2_data_2rxs   = join( in_dir, '2RXS_AllWISE_catalog_paper_2017May26_v5_11_0_sdss_26_VERON_MASKED.fits'   )

data_file = path_2_data_2rxs

path_2_randoms = '/data44s/eroAGN_WG_DATA/DATA/randoms/randoms_2rxsmask_GAIA_star_mask.fits'

suffix = 'mask_gaia_g_lt_6'

# sub catalogs
# - full cat
# - cut abs(gal_lat)>20 degrees
# - cut on X-ray depth map
# - cut on wise depth map

# write ascii catalogs + param files here :
# 

ra_name_data = 'ALLW_RA'
dec_name_data = 'ALLW_DEC'

ra_name_rds = 'RA'
dec_name_rds = 'DEC'

hduD     = fits.open(data_file)
ra_data    = hduD[1].data[ra_name_data]
dec_data    = hduD[1].data[dec_name_data]
z_data = np.zeros_like(ra_data)
ratelim_data    = hduD[1].data['RATELIM']

coords = SkyCoord(ra_data, dec_data, unit='deg', frame='icrs')
bb_data = coords.galactic.b.value
ll_data = coords.galactic.l.value
bb_ecl_data = coords.barycentrictrueecliptic.lat

near_a_star_D = (hduD[1].data['nearby_gaia_star_-10_g_5'] |
	hduD[1].data['nearby_gaia_star_5_g_6']   |
	#hduD[1].data['nearby_gaia_star_6_g_7']   |
	#hduD[1].data['nearby_gaia_star_7_g_8']   |
	#hduD[1].data['nearby_gaia_star_8_g_9']   |
	#hduD[1].data['nearby_gaia_star_9_g_10']  |
	#hduD[1].data['nearby_gaia_star_10_g_11'] |
	#hduD[1].data['nearby_gaia_star_11_g_12'] |
	#hduD[1].data['nearby_gaia_star_12_g_13'] |
	#hduD[1].data['nearby_gaia_star_13_g_14'] |
	#hduD[1].data['nearby_gaia_star_14_g_15'] |
	#hduD[1].data['nearby_gaia_star_15_g_16'] |
	#hduD[1].data['nearby_gaia_star_16_g_17'] 
	)


stars_data = (hduD[1].data['p_any']>0.5)&(hduD[1].data['2RXS_ExiML']>10)
rt_sel_data = (ratelim_data>0.015)
x_gal_data = (abs(bb_data)>20)&(dec_data<80)&(dec_data>-80)&(bb_ecl_data.value>-80)
selection_data = (x_gal_data)&(stars_data==False)&(rt_sel_data) &(near_a_star_D==False)

N_data = len(ra_data[selection_data])

#print( hduD[1].data.dtype)

hduR     = fits.open(path_2_randoms)
ra_rds    = hduR[1].data[ra_name_rds]
dec_rds    = hduR[1].data[dec_name_rds]
ratelim_rds    = hduR[1].data['RATELIM']

rt_sel_data = (ratelim_rds>0.015)

coords = SkyCoord(ra_rds, dec_rds, unit='deg', frame='icrs')
bb_rds = coords.galactic.b.value
ll_rds = coords.galactic.l.value
bb_ecl_rds = coords.barycentrictrueecliptic.lat

x_gal_rds = (abs(bb_rds)>20)&(dec_rds<80)&(dec_rds>-80)&(bb_ecl_rds.value>-80)

N_rds = len(ra_rds[x_gal_rds])

rd_rds = np.random.random(len(ra_rds))

# random downsampling to obtain 20x the randoms
frac = 20./(N_rds/N_data)
down_samp = (rd_rds<frac)


near_a_star_R = (hduR[1].data['nearby_gaia_star_-10_g_5'] |
	hduR[1].data['nearby_gaia_star_5_g_6']   |
	#hduR[1].data['nearby_gaia_star_6_g_7']   |
	#hduR[1].data['nearby_gaia_star_7_g_8']   |
	#hduR[1].data['nearby_gaia_star_8_g_9']   |
	#hduR[1].data['nearby_gaia_star_9_g_10']  |
	#hduR[1].data['nearby_gaia_star_10_g_11'] |
	#hduR[1].data['nearby_gaia_star_11_g_12'] |
	#hduR[1].data['nearby_gaia_star_12_g_13'] |
	#hduR[1].data['nearby_gaia_star_13_g_14'] |
	#hduR[1].data['nearby_gaia_star_14_g_15'] |
	#hduR[1].data['nearby_gaia_star_15_g_16'] |
	#hduR[1].data['nearby_gaia_star_16_g_17'] 
	)


selection_rds = (x_gal_rds) & (down_samp) & (rt_sel_data)  &(near_a_star_R==False)

N_rds = len(ra_rds[selection_rds])

out_data = os.path.join(out_dir , '2RXS_AllWISE_catalog_paper_2017May26_X_GAL_'+suffix+'.data')
np.savetxt(out_data, np.transpose([ra_data[selection_data], dec_data[selection_data], 0.7*np.ones_like(ra_data[selection_data]), np.ones_like(ra_data[selection_data]) ])  )
print(out_data)

out_rds = os.path.join(out_dir , '2RXS_AllWISE_catalog_paper_2017May26_X_GAL_'+suffix+'.random')
np.savetxt(out_rds, np.transpose([ra_rds[selection_rds], dec_rds[selection_rds], 0.7*np.ones_like(ra_rds[selection_rds]), np.ones_like(ra_rds[selection_rds]) ])  )
print(out_rds)

