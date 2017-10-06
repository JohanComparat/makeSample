import astropy.io.fits as fits
import os
import sys
from os.path import join
import pymangle as mangle
import numpy as np

def mask_catalog(in_file, out_file, ra_name, dec_name):
	print in_file
	hdu     = fits.open(in_file)
	cols    = hdu[1].columns
	data    = hdu[1].data
	ra      = data[ra_name]
	dec     = data[dec_name]

	# masks
	maskdir    = join(os.environ['OBS_REPO'],'SDSS/dr14/spiders/masks')
	dict_mask  = {}
	dict_mask['mask_bright_object_rykoff'] = join(maskdir,'bright_object_mask_rykoff_pix.ply')
	dict_mask['mask_Tycho20Vmag10']   = join(maskdir, 'tycho2mask-0Vmag10.pol')
	dict_mask['mask_Tycho210Vmag11']  = join(maskdir, 'tycho2mask-10Vmag11.pol')
	dict_mask['mask_Tycho211Vmag115'] = join(maskdir, 'tycho2mask-11Vmag115.pol')
	dict_mask['mask_dr14_N'] = join(maskdir, 'mask-QSO-N-eboss_v1.9f4.ply')
	dict_mask['mask_dr14_S'] = join(maskdir, 'mask-QSO-S-eboss_v1.9f4.ply')
	dict_mask['unphot-ugriz'] = join(maskdir, 'badfield_mask_unphot-ugriz_pix.ply')

	def createBooleanColumn(name, ra, dec, dict_mask=dict_mask):
		#print dict_mask[name]
		mng         = mangle.Mangle(dict_mask[name])
		polyid, w   = mng.polyid_and_weight(ra,dec)
		#print polyid, w
		tmp         = (polyid!=-1)
		col = fits.Column(name=name, format="L", array=tmp)
		colw = fits.Column(name=name+"_w", format="D", array=w)
		return col, colw

	# Create the mangle mask columns: bright star mask and bright object mask :
	for mask in dict_mask.keys():
		print mask#, ra, dec, dict_mask
		c1, c2 = createBooleanColumn(mask, ra, dec)
		cols += c1
		cols += c2
		

	# writes 
	hdu2 = fits.BinTableHDU.from_columns(cols)
	if os.path.isfile(out_file):
	   os.system('rm '+out_file)
	   
	hdu2.writeto(out_file)



in_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask.fits')
out_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask_DR14area.fits')
mask_catalog(in_file, out_file, 'ALLW_RA', 'ALLW_DEC')

in_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/randoms/randoms_dr14_NGC_with2RXS_mask.fits')
out_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/randoms/randoms_dr14_NGC_with2RXS_mask_DR14area.fits')
mask_catalog(in_file, out_file, 'RA', 'DEC')

in_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/randoms/randoms_dr14_SGC_with2RXS_mask.fits')
out_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/randoms/randoms_dr14_SGC_with2RXS_mask_DR14area.fits')
mask_catalog(in_file, out_file, 'RA', 'DEC')

