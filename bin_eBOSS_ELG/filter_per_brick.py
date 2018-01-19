"""
Contributions by A. Raichoor
"""
import numpy as n
import astropy.io.fits as fits
import os
import sys
from os.path import join

# survey bricks
hdu_br=fits.open(join(os.environ['OBS_REPO'], 'SDSS', 'targets','decals','survey-bricks.fits.gz'))
data=hdu_br[1].data
for quant in ['BRICKNAME','RA1', 'RA2','DEC1','DEC2']:
    exec quant+' = data["'+quant+'"]'

# dr5-eboss2 brick list
data_file_dr5_2  = join(os.environ['OBS_REPO'], 'SDSS', 'targets', 'ELG', 'elg_190_ngc_DR5_2.mask.eboss25.fits')
hdu=fits.open(data_file_dr5_2)
brn_uniq = n.unique(hdu[1].data['brickname'])

#data1/SDSS/targets/ELG/random_190_ngc_DR5_2.mask.fits
#data1/SDSS/targets/ELG/random_190_ngc_DR5.mask.fits

# randoms
randfits  = join(os.environ['OBS_REPO'], 'SDSS', 'targets', 'ELG', 'random_190_ngc_DR5.mask.fits')
hdu_R = fits.open(randfits)
ra  = hdu_R[1].data['ra']
dec = hdu_R[1].data['dec']
keep= n.zeros(len(ra),dtype='bool')
for brn in brn_uniq:
    tmp = (BRICKNAME==brn)
    ra1 = RA1[tmp][0]
    ra2 = RA2[tmp][0]
    dec1= DEC1[tmp][0]
    dec2= DEC2[tmp][0]
    tmp = (ra>ra1) & (ra<ra2) & (dec>dec1) & (dec<dec2)
    keep[tmp] = True
    
out_file = join(os.environ['OBS_REPO'], 'SDSS', 'targets', 'ELG', 'random_190_ngc_DR5.mask.common_bricks.fits')
new_columns = hdu_R[1].data.columns 
hdu2 = fits.BinTableHDU.from_columns(new_columns)
hdu2.data = hdu2.data[keep]
if os.path.isfile(out_file):
	os.system("rm "+out_file)

hdu2.writeto(out_file ) # for DR2 
