import astropy.io.fits as fits
import os
import sys
from os.path import join
#import pymangle as mangle
import numpy as np
import matplotlib.pyplot as p
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

import astropy.units as u
import astropy.cosmology as cc
cosmo = cc.Planck13

data_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask_DR14area_SPM.fits')
rds_n_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/randoms/randoms_dr14_NGC_with2RXS_mask_DR14area.fits')
rds_s_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/randoms/randoms_dr14_SGC_with2RXS_mask_DR14area.fits')

rds_file = rds_n_file

ra_name_data = 'ALLW_RA'
dec_name_data = 'ALLW_DEC'

ra_name_rds = 'RA'
dec_name_rds = 'DEC'

z_max = 1.3

#def run_ks_test(data_file, ra_name_data, dec_name_data, rds_file, ra_name_rds, dec_name_rds, ratelim_min):

hduD     = fits.open(data_file)
ra_data    = hduD[1].data[ra_name_data]
dec_data    = hduD[1].data[dec_name_data]
z_data = np.zeros_like(ra_data)
ratelim_data    = hduD[1].data['MASK_2RXS_RATELIM']

hduR     = fits.open(rds_file)
ra_rds    = hduR[1].data[ra_name_rds]
dec_rds    = hduR[1].data[dec_name_rds]
ratelim_rds    = hduR[1].data['MASK_2RXS_RATELIM']
down_samp = (np.random.random(len(ra_rds))<0.1)

def plot_clusterings_catalog(out_file, ratelim = 0.0, hemisphere='N'):
  spec_SDSS = (hduD[1].data['DR14_Z']>0)&(hduD[1].data['DR14_Z_ERR']>0)&(hduD[1].data['DR14_ZWARNING']==0)
  spec_veron = (hduD[1].data['z_veron']>0)
  z_data[spec_veron] = hduD[1].data['z_veron'][spec_veron]
  z_data[spec_SDSS] = hduD[1].data['DR14_Z'][spec_SDSS]
  spectro = (z_data>0)
  starsD = (hduD[1].data['mask_Tycho20Vmag10']==False)& (hduD[1].data['mask_Tycho210Vmag11']==False)& (hduD[1].data['mask_bright_object_rykoff']==False) # (hduD[1].data['unphot-ugriz']==False)& 
  starsR = (hduR[1].data['mask_Tycho20Vmag10']==False) &(hduR[1].data['mask_Tycho210Vmag11']==False) &(hduR[1].data['mask_bright_object_rykoff']==False) # (hduR[1].data['unphot-ugriz']==False)& 

  if hemisphere=='N':
    sel_data = (z_data<z_max)&(ratelim_data>ratelim) &(hduD[1].data['mask_dr14_N'])  &(hduD[1].data['mask_dr14_N_w']>0.9) &(spectro)&(starsD) 
    sel_rds = (down_samp)&(ratelim_rds>ratelim) &(hduR[1].data['mask_dr14_N']) &(hduR[1].data['mask_dr14_N_w']>0.9) &(starsR) 
  if hemisphere=='S':
    sel_data = (z_data<z_max)&(ratelim_data>ratelim) &(hduD[1].data['mask_dr14_S'])  &(hduD[1].data['mask_dr14_S_w']>0.9) &(spectro)&(starsD) 
    sel_rds = (down_samp)&(ratelim_rds>ratelim) &(hduR[1].data['mask_dr14_S']) &(hduR[1].data['mask_dr14_S_w']>0.9) &(starsR)
  
  N_data = len(ra_data[sel_data]) 
  N_rds = len(ra_rds[sel_rds]) 
  print(N_data, N_rds)
  dz=0.05
  zs=np.arange(0.,1.4, dz)
  nn,bb = np.histogram(z_data[sel_data], bins=zs)#, weights=1./w_col.array)
  nz=interp1d((zs[1:]+zs[:-1])/2.,nn)
  rdsz=[]
  for i in range(1,len(zs)-1,1):
    inter=np.random.uniform(low=zs[i]-dz/2., high=zs[i]+dz/2., size=int( 1000* nz( zs[i] )))
    rdsz.append(inter)

  rds=np.hstack((rdsz))
  np.random.shuffle(rds)
  print(len(rds))
  RR=rds[:N_rds]#-dz/2.

  p.figure(1, (5,5))
  p.axes([0.17,0.17,0.78,0.78])
  p.plot(ra_data[sel_data], dec_data[sel_data], 'k+')
  p.xlabel('RA [deg]')
  p.ylabel('DEC [deg]')
  p.title('data dr14')
  p.savefig(out_file+str(ratelim)+"_data_sky.png")
  p.clf()

  p.figure(1, (5,5))
  p.axes([0.17,0.17,0.78,0.78])
  p.plot(ra_rds[sel_rds], dec_rds[sel_rds], 'k,')
  p.xlabel('RA [deg]')
  p.ylabel('DEC [deg]')
  p.title('random dr14')
  p.savefig(out_file+str(ratelim)+"_randoms_sky.png")
  p.clf()

  flux = hduD[1].data['2RXS_SRC_FLUX'][sel_data]

  p.figure(1, (5,5))
  p.axes([0.17,0.17,0.78,0.78])
  p.plot(z_data[sel_data], hduD[1].data['2RXS_SRC_FLUX'][sel_data], 'k,')
  p.xlabel('redshift')
  p.ylabel('Xray flux [0.5-2 keV, erg/cm2/s]')
  p.yscale('log')
  p.ylim((1e-13, 5e-11))
  p.title('data dr14')
  p.savefig(out_file+str(ratelim)+"_data_Xray_flux_redshift.png")
  p.clf()

  dl_cm = cosmo.luminosity_distance(z_data[sel_data]).to(u.cm)
  LX = dl_cm.value**2. * flux * 4 * np.pi
  LX_125 = dl_cm.value**2.* 10**(-12.5)* 4 * np.pi
  LX_127 = dl_cm.value**2.* 10**(-12.7)* 4 * np.pi

  p.figure(1, (5,5))
  p.axes([0.17,0.17,0.78,0.78])
  p.plot(z_data[sel_data], LX, 'k+', label='data')
  p.plot(z_data[sel_data], LX_125, 'r,', label=r'$10^{-12.5}$erg/cm2/s')
  p.plot(z_data[sel_data], LX_127, 'b,', label=r'$10^{-12.7}$erg/cm2/s')
  p.xlabel('redshift')
  p.ylabel('Xray luminosity [0.5-2 keV, erg/s] ')
  p.yscale('log')
  p.ylim((1e41, 1e47))
  p.legend(frameon=False, loc=0)
  p.title('data dr14')
  p.savefig(out_file+str(ratelim)+"_data_Xray_luminosity_redshift.png")
  p.clf()


print(' / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / ')
print('NORTH')
out_file  = join(os.environ['HOME'], 'wwwDir/eRoMok/clustering/data/clustering_agn_N_RL_')

plot_clusterings_catalog(out_file+'_N', ratelim = 0.005, hemisphere='N')
plot_clusterings_catalog(out_file+'_N', ratelim = 0.015, hemisphere='N')

rds_file = rds_s_file

hduR     = fits.open(rds_file)
ra_rds    = hduR[1].data[ra_name_rds]
dec_rds    = hduR[1].data[dec_name_rds]
ratelim_rds    = hduR[1].data['MASK_2RXS_RATELIM']
down_samp = (np.random.random(len(ra_rds))<0.1)

plot_clusterings_catalog(out_file+'_S', ratelim = 0.005, hemisphere='S')
plot_clusterings_catalog(out_file+'_S', ratelim = 0.015, hemisphere='S')


sys.exit()

print(' / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / ')
print('SOUTH')
rds_file = rds_s_file

hduR     = fits.open(rds_file)
ra_rds    = hduR[1].data[ra_name_rds]
dec_rds    = hduR[1].data[dec_name_rds]
ratelim_rds    = hduR[1].data['MASK_2RXS_RATELIM']
down_samp = (np.random.random(len(ra_rds))<0.1)

