"""
removes from dr5-eboss all bricks present in dr5-eboss2
"""
import numpy as n
import astropy.io.fits as fits
import os
import sys
from os.path import join

data_file_dr5_2  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/elg_190_ngc_DR5_2.mask.eboss25.fits')
data_file_dr5  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/elg_190_ngc_DR5.mask.eboss25.fits')
out_file = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/elg_190_ngc_DR5.mask.eboss25.nocommonbricks.fits')

d52     = fits.open(data_file_dr5_2)[1].data
d5     = fits.open(data_file_dr5)[1].data

b52 = n.array(list(set(d52['brickname'])))
b5 = n.array(list(set(d5['brickname'])))

common_bricks = n.in1d(d5['brickname'], b52)

sel_data = (common_bricks==False) 	
new_columns = d5.columns 
hdu2 = fits.BinTableHDU.from_columns(new_columns)
hdu2.data = hdu2.data[sel_data]
if os.path.isfile(out_file):
	os.system("rm "+out_file)

hdu2.writeto(out_file ) # for DR2 

"""
out_file = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/elg_190_ngc_DR5.mask.eboss25.samebricks.fits')

d52     = fits.open(data_file_dr5_2)[1].data
d5     = fits.open(data_file_dr5)[1].data

b52 = n.array(list(set(d52['brickname'])))
b5 = n.array(list(set(d5['brickname'])))

common_bricks = n.in1d(d5['brickname'], b52)

sel_data = (common_bricks==True) 	
new_columns = d5.columns 
hdu2 = fits.BinTableHDU.from_columns(new_columns)
hdu2.data = hdu2.data[sel_data]
if os.path.isfile(out_file):
	os.system("rm "+out_file)

hdu2.writeto(out_file ) # for DR2 
"""