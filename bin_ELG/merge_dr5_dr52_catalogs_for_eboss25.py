import numpy as n
import astropy.io.fits as fits
import os
import sys
from os.path import join

def write_target_cat(data_file, out_file):
	ra_name_data = 'ra'
	dec_name_data = 'dec'
	dat     = fits.open(data_file)[1].data
	ra_data    = dat[ra_name_data]
	dec_data    = dat[dec_name_data]
	starsD = (dat['mask_Tycho20Vmag10']==False)& (dat['mask_Tycho210Vmag11']==False)& (dat['mask_bright_object_rykoff']==False) & (dat['mask_Tycho211Vmag115']==False) & (dat['mask_Tycho2115Vmag120']==False)
	
	sel_data = (starsD) & (dat['eboss25']) 	
	new_columns = dat.columns 
	hdu2 = fits.BinTableHDU.from_columns(new_columns)
	hdu2.data = hdu2.data[sel_data]
	if os.path.isfile(out_file):
		os.system("rm "+out_file)

	hdu2.writeto(out_file ) # for DR2 


data_file_dr5_2  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/elg_190_ngc_DR5_2.mask.eboss25.fits')
data_file_dr5  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/elg_190_ngc_DR5.mask.eboss25.fits')

	d52     = fits.open(data_file_dr5_2)[1].data
	d5     = fits.open(data_file_dr5)[1].data

	b52 = n.array(list(set(d52['brickname'])))
	b5 = n.array(list(set(d5['brickname'])))



write_target_cat(data_file, data_file[:-5]+".eboss25.fits")

write_target_cat(data_file, data_file[:-5]+".eboss25.fits")

