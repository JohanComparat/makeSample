import astropy.io.fits as fits
import os
import sys
from os.path import join
#import pymangle as mangle
import numpy as np
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
data_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask_DR14area.fits')
rds_n_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/randoms/randoms_dr14_NGC_with2RXS_mask_DR14area.fits')
rds_s_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/randoms/randoms_dr14_SGC_with2RXS_mask_DR14area.fits')

rds_file = rds_n_file

ra_name_data = 'ALLW_RA'
dec_name_data = 'ALLW_DEC'

ra_name_rds = 'RA'
dec_name_rds = 'DEC'

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

def write_ascii_clusterings_catalog(out_file, ratelim_min = 0.01):
  spec_SDSS = (hduD[1].data['DR14_Z']>0)&(hduD[1].data['DR14_Z_ERR']>0)&(hduD[1].data['DR14_ZWARNING']==0)
  spec_veron = (hduD[1].data['z_veron']>0)
  z_data[spec_veron] = hduD[1].data['z_veron'][spec_veron]
  z_data[spec_SDSS] = hduD[1].data['DR14_Z'][spec_SDSS]
  spectro = (z_data>0)
  starsD = (hduD[1].data['mask_Tycho20Vmag10']==False)& (hduD[1].data['mask_Tycho210Vmag11']==False)& (hduD[1].data['mask_bright_object_rykoff']==False) # (hduD[1].data['unphot-ugriz']==False)& 
  sel_data = (ratelim_data>ratelim_min) &(hduD[1].data['mask_dr14_N'])  &(hduD[1].data['mask_dr14_N_w']>0.9) &(spectro)&(starsD) #&(ra_data>140.)&(ra_data<300.)
  starsR = (hduR[1].data['mask_Tycho20Vmag10']==False) &(hduR[1].data['mask_Tycho210Vmag11']==False) &(hduR[1].data['mask_bright_object_rykoff']==False) # (hduR[1].data['unphot-ugriz']==False)& 
  sel_rds = (down_samp)&(ratelim_rds>ratelim_min) &(hduR[1].data['mask_dr14_N']) &(hduR[1].data['mask_dr14_N_w']>0.9) &(starsR) #&(ra_rds>140.)&(ra_rds<300.)
  # ra figure
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

  np.savetxt(out_file+str(ratelim_min)+".data", np.transpose([ ra_data[sel_data], dec_data[sel_data], z_data[sel_data], np.ones_like(z_data[sel_data]) ]))
  np.savetxt(out_file+str(ratelim_min)+".random", np.transpose([ ra_rds[sel_rds], dec_rds[sel_rds], RR, np.ones_like(RR) ]))

print(' / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / ')
print('NORTH')
out_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target/clustering_agn_N_RL_')

write_ascii_clusterings_catalog(out_file, ratelim_min = 0.005)
write_ascii_clusterings_catalog(out_file, ratelim_min = 0.015)
write_ascii_clusterings_catalog(out_file, ratelim_min = 0.0216)

#sys.exit()

print(' / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / ')
print('SOUTH')
rds_file = rds_s_file

hduR     = fits.open(rds_file)
ra_rds    = hduR[1].data[ra_name_rds]
dec_rds    = hduR[1].data[dec_name_rds]
ratelim_rds    = hduR[1].data['MASK_2RXS_RATELIM']
down_samp = (np.random.random(len(ra_rds))<0.1)

def write_ascii_clusterings_catalog(out_file, ratelim_min = 0.01):
  spec_SDSS = (hduD[1].data['DR14_Z']>0)&(hduD[1].data['DR14_Z_ERR']>0)&(hduD[1].data['DR14_ZWARNING']==0)
  spec_veron = (hduD[1].data['z_veron']>0)
  z_data[spec_veron] = hduD[1].data['z_veron'][spec_veron]
  z_data[spec_SDSS] = hduD[1].data['DR14_Z'][spec_SDSS]
  spectro = (z_data>0)
  starsD = (hduD[1].data['mask_Tycho20Vmag10']==False)& (hduD[1].data['mask_Tycho210Vmag11']==False)& (hduD[1].data['mask_bright_object_rykoff']==False) # (hduD[1].data['unphot-ugriz']==False)& 
  sel_data = (ratelim_data>ratelim_min) &(hduD[1].data['mask_dr14_S'])  &(hduD[1].data['mask_dr14_S_w']>0.9) &(spectro)&(starsD) #&(ra_data>140.)&(ra_data<300.)
  starsR = (hduR[1].data['mask_Tycho20Vmag10']==False) &(hduR[1].data['mask_Tycho210Vmag11']==False) &(hduR[1].data['mask_bright_object_rykoff']==False) # (hduR[1].data['unphot-ugriz']==False)& 
  sel_rds = (down_samp)&(ratelim_rds>ratelim_min) &(hduR[1].data['mask_dr14_S']) &(hduR[1].data['mask_dr14_S_w']>0.9) &(starsR) #&(ra_rds>140.)&(ra_rds<300.)
  # ra figure
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

  np.savetxt(out_file+str(ratelim_min)+".data", np.transpose([ ra_data[sel_data], dec_data[sel_data], z_data[sel_data], np.ones_like(z_data[sel_data]) ]))
  np.savetxt(out_file+str(ratelim_min)+".random", np.transpose([ ra_rds[sel_rds], dec_rds[sel_rds], RR, np.ones_like(RR) ]))

  
out_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target/clustering_agn_S_RL_')

write_ascii_clusterings_catalog(out_file, ratelim_min = 0.005)
write_ascii_clusterings_catalog(out_file, ratelim_min = 0.015)
write_ascii_clusterings_catalog(out_file, ratelim_min = 0.0216)
