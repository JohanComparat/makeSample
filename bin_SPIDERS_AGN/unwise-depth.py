# compute unwise depth


import healpy
import numpy as np
import os, sys
import astropy.io.fits as fits


nl = lambda selection : len(selection.nonzero()[0])
# unWISE file list
top_dir = '/data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/unwise/release/'
file_list = np.loadtxt(os.path.join(top_dir, 'cat.list'), unpack=True, dtype='str')

nside_values = 2**np.arange(10,11,1)

for ii in range(len(file_list)):
	out_file = os.path.join(top_dir, 'depth-mask', 'mask-'+file_list[ii])
	f0 = fits.open(os.path.join(top_dir, 'band-merged', file_list[ii]))

	N_targets_total = len(f0[1].data['ra'])


	HEALPIX_VAL={}
	set_HP = {}
	area_estimated = {}
	density_estimated = {}
	number_per_pixel = {}
	density_per_pixel = {}
	median_dfluxlbs_per_pixel = {}
	median_fwhm_per_pixel = {}
	ra_hp = {}
	dec_hp = {}

	prihdr = fits.Header()
	prihdr['author'] = 'COMPARAT'
	prihdu = fits.PrimaryHDU(header=prihdr)

	hdu_list = [prihdu]
	for nside in nside_values:
		HEALPIX_VAL[nside] = healpy.ang2pix(nside,  f0[1].data['dec']*np.pi/180.+ np.pi/2. , f0[1].data['ra']*np.pi/180. , nest=True)
		set_HP[nside] = np.array(list(set(HEALPIX_VAL[nside]   )))
		area_estimated[nside] = len( set_HP[nside]    ) *healpy.nside2pixarea(nside  , degrees=True)
		number_per_pixel[nside]    = np.array([nl((HEALPIX_VAL[nside] == el)) for el in set_HP[nside]])  
		density_estimated[nside]    = N_targets_total / area_estimated[nside]   
		density_per_pixel[nside]    = number_per_pixel[nside]    /healpy.nside2pixarea(nside , degrees=True)
		median_dfluxlbs_per_pixel[nside]    = np.array([np.median(f0[1].data['dfluxlbs'][(HEALPIX_VAL[nside] == el)]) for el in set_HP[nside]])  
		median_fwhm_per_pixel[nside]    = np.array([np.median(f0[1].data['fwhm'][(HEALPIX_VAL[nside] == el)]) for el in set_HP[nside]])  
		ra_hp = np.array([ healpy.pix2ang(nside, pix_id, nest=True)[1]*180./np.pi for pix_id in set_HP[nside] ])
		dec_hp = np.array([ (healpy.pix2ang(nside, pix_id, nest=True)[0]-np.pi/2.)*180./np.pi for pix_id in set_HP[nside] ])

	for nside in nside_values:
		hdu_cols = fits.ColDefs([
		fits.Column(name='HP_IDX'            , format='K' ,  array = HEALPIX_VAL[nside] ),   
		fits.Column(name='N'                 , format='K' , array = number_per_pixel[nside] ),   
		fits.Column(name='RHO'               , format='1E' ,  array = density_per_pixel[nside] ),   
		fits.Column(name='MED_DFLUXLBS'      , format='1E' ,  array = median_dfluxlbs_per_pixel[nside] ),   
		fits.Column(name='MED_FWHM'          , format='1E' ,  array = median_fwhm_per_pixel[nside] ),  
		fits.Column(name='RA_PIX'          , format='1E' ,  array = ra_hp[nside] ),  
		fits.Column(name='DEC_PIX'          , format='1E' ,  array = dec_hp[nside] )
		])
		hdu = fits.BinTableHDU.from_columns(hdu_cols)
		hdu.name = 'healpix_'+str(nside)
		hdu.header['NSIDE'] = nside
		hdu_list.append( hdu )


	#writes the file
	thdulist = fits.HDUList(hdu_list)

	print( out_file )
	if os.path.isfile(out_file):
		os.system("rm "+out_file)
	thdulist.writeto(out_file)



#for nside in nside_values:
	#print('=======================', nside,'==============================')
	#print(number_per_pixel[nside].min(), number_per_pixel[nside].max(), number_per_pixel[nside].mean())
	#print(density_per_pixel[nside].min(), density_per_pixel[nside].max(), density_per_pixel[nside].mean())
	#print(median_fwhm_per_pixel[nside].min(), median_fwhm_per_pixel[nside].max(), median_fwhm_per_pixel[nside].mean())
	#print(median_dfluxlbs_per_pixel[nside].min(), median_dfluxlbs_per_pixel[nside].max(), median_dfluxlbs_per_pixel[nside].mean())

#top_dir = '/data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/unwise/release/'
#file_list = np.loadtxt(os.path.join(top_dir, 'cat.list'), unpack=True, dtype='str')

#top_dir = /data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/unwise/release/band-merged

#add healpix 
#take medians in pixels for 


#dflux
#dfluxlbs
#fwhm


#(0)comparat@ds52:~/eroAGN_WG_DATA/DATA$ fstruct photometry/catalogs/unwise/release/band-merged/3584p196.cat.fits
  #No. Type     EXTNAME      BITPIX Dimensions(columns)      PCOUNT  GCOUNT
 
   #0  PRIMARY                  8     0                           0    1
   #1  BINTABLE                 8     262(27) 82512               0    1
 
      #Column Name                Format     Dims       Units     TLMIN  TLMAX
      #1 x                          2D
      #2 y                          2D
      #3 flux                       2E
      #4 dx                         2E
      #5 dy                         2E
      #6 dflux                      2E
      #7 qf                         2E
      #8 rchi2                      2E
      #9 fracflux                   2E
     #10 fluxlbs                    2E
     #11 dfluxlbs                   2E
     #12 fwhm                       2E
     #13 spread_model               2E
     #14 dspread_model              2E
     #15 sky                        2E
     #16 ra12                       2D
     #17 dec12                      2D
     #18 coadd_id                   8A
     #19 unwise_detid               36A     (18,2)
     #20 nm                         2I
     #21 primary12                  2I
     #22 flags_unwise               2I
     #23 flags_info                 2I
     #24 ra                         D
     #25 dec                        D
     #26 primary                    I
     #27 unwise_objid               16A
 
