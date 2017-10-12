import astropy.io.fits as fits
import os
import sys
from os.path import join
import pymangle as mangle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p


data_file  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/elg_190_ngc_DR5.mask.fits')
rds_file  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/random_190_ngc_DR5.mask.fits')

data_file  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/elg_190_ngc_DR5_2.mask.fits')
rds_file  = join(os.environ['OBS_REPO'], 'SDSS', 'targets/ELG/random_190_ngc_DR5_2.mask.fits')


ra_name_data = 'ra'
dec_name_data = 'dec'

ra_name_rds = 'ra'
dec_name_rds = 'dec'

#def run_ks_test(data_file, ra_name_data, dec_name_data, rds_file, ra_name_rds, dec_name_rds, ratelim_min):

hduD     = fits.open(data_file)
ra_data    = hduD[1].data[ra_name_data]
dec_data    = hduD[1].data[dec_name_data]
z_data = np.zeros_like(ra_data)

hduR     = fits.open(rds_file)
ra_rds    = hduR[1].data[ra_name_rds]
dec_rds    = hduR[1].data[dec_name_rds]

areaD = (ra_data>135.)&(ra_data<170.)&(dec_data>21.5)&(dec_data<33)
starsD = (hduD[1].data['mask_Tycho20Vmag10']==False)& (hduD[1].data['mask_Tycho210Vmag11']==False)& (hduD[1].data['mask_bright_object_rykoff']==False) & (hduD[1].data['mask_Tycho211Vmag115']==False) 
sel_data = (areaD) & (starsD) 

areaR = (ra_rds>135.)&(ra_rds<170.)&(dec_rds>21.5)&(dec_rds<33)
starsR = (hduR[1].data['mask_Tycho20Vmag10']==False) & (hduR[1].data['mask_Tycho210Vmag11']==False) & (hduR[1].data['mask_bright_object_rykoff']==False) & (hduR[1].data['mask_Tycho211Vmag115']==False)

sel_rds = (areaR)&(starsR)

# ra figure
p.figure(1, (5,5))
N_ra_data, bb, ptxs = p.hist(ra_data[sel_data], normed=True, cumulative=True, bins=1000, histtype='step', label='data')
N_ra_rds, bb, ptxs = p.hist(ra_rds[sel_rds], normed=True, cumulative=True, bins=bb, histtype='step', label='rds')
p.xlabel('ra')
p.ylabel('normed cumulative distribution')
p.ylim((-0.01,1.01))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig(join(os.environ['HOME'], 'wwwDir/sdss/elg/N_ra_ks_test.png'))
p.clf()
# dec figure
p.figure(1, (5,5))
N_dec_data, bb, ptxs = p.hist(dec_data[sel_data], normed=True, cumulative=True, bins=1000, histtype='step', label='data ')
N_dec_rds, bb, ptxs = p.hist(dec_rds[sel_rds], normed=True, cumulative=True, bins=bb, histtype='step', label='rds ')
p.xlabel('dec')
p.ylabel('normed cumulative distribution')
p.ylim((-0.01,1.01))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig(join(os.environ['HOME'], 'wwwDir/sdss/elg/N_dec_ks_test.png'))
p.clf()
# ra-dec figure
p.figure(1, (5,5))
p.plot(ra_rds[sel_rds], dec_rds[sel_rds], 'k,', rasterized=True)
p.xlabel('ra')
p.ylabel('dec')
#p.legend(loc=0, frameon=False)
p.savefig(join(os.environ['HOME'], 'wwwDir/sdss/elg/N_sky_ks_test_D.png'))
p.clf()

p.figure(1, (5,5))
p.plot(ra_data[sel_data], dec_data[sel_data], 'r,', rasterized=True)
p.xlabel('ra')
p.ylabel('dec')
#p.legend(loc=0, frameon=False)
p.savefig(join(os.environ['HOME'], 'wwwDir/sdss/elg/N_sky_ks_test_R.png'))
p.clf()

print('=====================================================================')
N_data = len(ra_data[sel_data]) 
N_rds = len(ra_rds[sel_rds]) 
print ( " N AGN", N_data    )
print ( " N randoms", N_rds )
Dnm_ra = np.max(N_ra_rds-N_ra_data)
Dnm_dec = np.max(N_dec_rds-N_dec_data)
# for the 0.05 level of rejection
Dlimit_value_005 = 1.36*((1.*N_data+1.*N_rds)/(1.*N_data*N_rds))**0.5
print( "Dnm limiting value", Dlimit_value_005)
print ( " KS-test ra" , Dnm_ra , "H rejected ?", Dnm_ra>Dlimit_value_005 )
print ( " KS-test dec", Dnm_dec, "H rejected ?", Dnm_dec>Dlimit_value_005 )




def write_ascii_clusterings_catalog(out_file_name, ratelim_min = 0.015):
  # first the data 
  spec_SDSS = (hduD[1].data['DR14_Z']>0)&(hduD[1].data['DR14_Z_ERR']>0)&(hduD[1].data['DR14_ZWARNING']==0)
  spec_veron = (hduD[1].data['z_veron']>0)
  z_data[spec_veron] = hduD[1].data['z_veron'][spec_veron]
  z_data[spec_SDSS] = duD[1].data['DR14_Z'][spec_SDSS]
  spectro = (z_data>0)
  starsD = (hduD[1].data['mask_Tycho20Vmag10']==False)& (hduD[1].data['mask_Tycho210Vmag11']==False)& (hduD[1].data['mask_bright_object_rykoff']==False) # (hduD[1].data['unphot-ugriz']==False)& 
  sel_data = (ratelim_data>ratelim_min) &(hduD[1].data['mask_dr14_N'])  &(hduD[1].data['mask_dr14_N_w']>0.9) &(spectro)&(starsD) #&(ra_data>140.)&(ra_data<300.)
  # what about the weights ?
  # correct from TSR ?
  # 
  np.savetxt(out_file_name, np.transpose([ra_data[sel_data], dec_data[sel_data], z_data[sel_data]]))
  # now the randoms 
  starsR = (hduR[1].data['mask_Tycho20Vmag10']==False) &(hduR[1].data['mask_Tycho210Vmag11']==False) &(hduR[1].data['mask_bright_object_rykoff']==False) # (hduR[1].data['unphot-ugriz']==False)& 
  sel_rds = (down_samp)&(ratelim_rds>ratelim_min) &(hduR[1].data['mask_dr14_N']) &(hduR[1].data['mask_dr14_N_w']>0.9) &(starsR) #&(ra_rds>140.)&(ra_rds<300.)
  #ra_rds[sel_rds], dec_rds[sel_rds],
  



# masks
maskdir    = join(os.environ['OBS_REPO'],'SDSS/dr14/spiders/masks')
dict_mask  = {}
dict_mask['mask_bright_object_rykoff'] = join(maskdir,'bright_object_mask_rykoff_pix.ply')
dict_mask['mask_Tycho20Vmag10']   = join(maskdir, 'tycho2mask-0Vmag10.pol')
dict_mask['mask_Tycho210Vmag11']  = join(maskdir, 'tycho2mask-10Vmag11.pol')
dict_mask['mask_Tycho211Vmag115'] = join(maskdir, 'tycho2mask-11Vmag115.pol')
dict_mask['mask_dr14_N'] = join(maskdir, 'mask-QSO-N-eboss_v1.9f4.ply')
dict_mask['mask_dr14_S'] = join(maskdir, 'mask-QSO-S-eboss_v1.9f4.ply')

def createBooleanColumn(name, ra, dec, dict_mask=dict_mask):
	#print dict_mask[name]
	mng         = mangle.Mangle(dict_mask[name])
	polyid      = mng.polyid(ra,dec)
	tmp         = (polyid!=-1)
	col = fits.Column(name=name, format="L", array=tmp)
	return col

# Create the mangle mask columns: bright star mask and bright object mask :
for mask in dict_mask.keys():
	print mask, ra, dec, dict_mask
	cols += createBooleanColumn(mask, ra, dec)

# writes 
hdu2 = fits.BinTableHDU.from_columns(cols)
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

