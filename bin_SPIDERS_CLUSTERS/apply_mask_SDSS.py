import astropy.io.fits as fits
import os
import sys
from os.path import join
import pymangle as mangle
import numpy as np

maskdir = '/home/comparat/data/CODEX'

dict_mask  = {}
dict_mask['mask_bad_field'] = join(maskdir, 'badfield_mask_unphot-ugriz_pix.ply')
dict_mask['mask_bright_star'] = join(maskdir, 'bright_star_mask_pix.ply')
#dict_mask['mask_Tycho20Vmag10']        = join(maskdir, 'tycho2mask-0Vmag10.pol'           )
#dict_mask['mask_Tycho210Vmag11']       = join(maskdir, 'tycho2mask-10Vmag11.pol'          )
#dict_mask['mask_Tycho211Vmag115']      = join(maskdir, 'tycho2mask-11Vmag115.pol'         )
#dict_mask['mask_Tycho2115Vmag120']     = join(maskdir, 'tycho2mask-115Vmag12.pol'         )

cluster_file = join(maskdir, 'cat_spiders.fits')
cluster_file_masked = join(maskdir, 'cat_spiders_masked_X.fits')

hdu     = fits.open(cluster_file)
cols    = hdu[1].columns
data    = hdu[1].data
ra      = data['RA']
dec     = data['DEC']
	
def createBooleanColumn(name, ra, dec, dict_mask=dict_mask):
	#print dict_mask[name]
	mng         = mangle.Mangle(dict_mask[name])
	polyid, w   = mng.polyid_and_weight(ra,dec)
	#print polyid, w
	tmp         = (polyid!=-1)
	col = fits.Column(name=name, format="L", array=tmp)
	colw = fits.Column(name=name+"_w", format="D", array=w)
	return col, colw

col_bad_field, col_bad_field_w = createBooleanColumn('mask_bad_field', ra, dec, dict_mask=dict_mask)
col_bright_star, col_bright_star_w = createBooleanColumn('mask_bright_star', ra, dec, dict_mask=dict_mask)

cols += col_bad_field
cols += col_bad_field_w
cols += col_bright_star
cols += col_bright_star_w

# writes 
hdu2 = fits.BinTableHDU.from_columns(cols)
if os.path.isfile(cluster_file_masked):
	os.system('rm '+cluster_file_masked)
	   
hdu2.writeto(cluster_file_masked)
