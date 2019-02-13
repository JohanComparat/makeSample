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

path_2_randoms = '/data44s/eroAGN_WG_DATA/DATA/randoms/randoms_2rxsmask.fits'

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

#stars_data = (hduD[1].data['p_any']>0.5)&(hduD[1].data['2RXS_ExiML']>10)
x_gal_data = (abs(bb_data)>20)&(dec_data<80)&(dec_data>-80)&(bb_ecl_data>-80)#&(stars_data==False)
selection_data = (x_gal_data)

N_data = len(ra_data[selection_data])

#print( hduD[1].data.dtype)

hduR     = fits.open(path_2_randoms)
ra_rds    = hduR[1].data[ra_name_rds]
dec_rds    = hduR[1].data[dec_name_rds]
ratelim_rds    = hduR[1].data['RATELIM']

coords = SkyCoord(ra_rds, dec_rds, unit='deg', frame='icrs')
bb_rds = coords.galactic.b.value
ll_rds = coords.galactic.l.value
bb_ecl_rds = coords.barycentrictrueecliptic.lat

x_gal_rds = (abs(bb_rds)>20)&(dec_rds<80)&(dec_rds>-80)&(bb_ecl_rds>-80)

N_rds = len(ra_rds[x_gal_rds])

rd_rds = np.random.random(len(ra_rds))

# random downsampling to obtain 20x the randoms
frac = 20./(N_rds/N_data)
down_samp = (rd_rds<frac)

selection_rds = (x_gal_rds) & (down_samp)

N_rds = len(ra_rds[selection_rds])


out_data = os.path.join(out_dir , '2RXS_AllWISE_catalog_paper_2017May26_X_GAL.data')
np.savetxt(out_data, np.transpose([ra_data[selection_data], dec_data[selection_data], 0.7*np.ones_like(ra_data[selection_data]), np.ones_like(ra_data[selection_data]) ])  )

out_rds = os.path.join(out_dir , '2RXS_AllWISE_catalog_paper_2017May26_X_GAL.random')
np.savetxt(out_rds, np.transpose([ra_rds[selection_rds], dec_rds[selection_rds], 0.7*np.ones_like(ra_rds[selection_rds]), np.ones_like(ra_rds[selection_rds]) ])  )


