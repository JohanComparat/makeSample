import astropy.io.fits as fits
import os
import sys
from os.path import join
import pymangle as mangle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p



def write_ascii_cat(data_file, out_file):
	ra_name_data = 'ra'
	dec_name_data = 'dec'
	hduD     = fits.open(data_file)
	ra_data    = hduD[1].data[ra_name_data]
	dec_data    = hduD[1].data[dec_name_data]
	z_data = np.zeros_like(ra_data)
	areaD = (ra_data>135.)&(ra_data<170.)&(dec_data>21.5)&(dec_data<33)
	starsD = (hduD[1].data['mask_Tycho20Vmag10']==False)& (hduD[1].data['mask_Tycho210Vmag11']==False)& (hduD[1].data['mask_bright_object_rykoff']==False) & (hduD[1].data['mask_Tycho211Vmag115']==False) 
	sel_data = (areaD) & (starsD) 
	np.savetxt(out_file, np.transpose([ra_data[sel_data], dec_data[sel_data],n.ones_like(dec_data[sel_data])*0.86,  n.ones_like(dec_data[sel_data]) ]))


data_file  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/elg_190_ngc_DR5.mask.fits')
write_ascii_cat(data_file, data_file[:-5]+"wtheta.ascii")
data_file  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/random_190_ngc_DR5.mask.fits')
write_ascii_cat(data_file, data_file[:-5]+"wtheta.ascii")
data_file  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/elg_190_ngc_DR5_2.mask.fits')
write_ascii_cat(data_file, data_file[:-5]+"wtheta.ascii")
data_file  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/random_190_ngc_DR5_2.mask.fits')
write_ascii_cat(data_file, data_file[:-5]+"wtheta.ascii")

