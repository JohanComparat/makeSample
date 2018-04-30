import astropy.io.fits as fits
import os
import sys
from os.path import join
#import pymangle as mangle
import numpy as np
import matplotlib.pyplot as p
from scipy.interpolate import interp1d

import astropy.units as u
import astropy.cosmology as cc
cosmo = cc.Planck13

data_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask_DR14area_SPM.fits')
rds_n_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/randoms/randoms_dr14_NGC_with2RXS_mask_DR14area.fits')
rds_s_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/randoms/randoms_dr14_SGC_with2RXS_mask_DR14area.fits')
#nR_NGC = 5135183. # random points
#nR_SGC = 3292159. # random points
#area = 2391. # unweighted
#rho_random = (nR_NGC+nR_SGC)/area # 4122.96 / deg2
#area_N = nR_NGC / rho_random
#area_S = nR_SGC / rho_random
#nR_NGC_written = 32015.
#nR_SGC_written = 21620.
#frac_NGC = 100.*nR_NGC_written/nR_NGC
#frac_SGC = 100.*nR_SGC_written/nR_SGC
#area_obs_N = frac_NGC * area_N # 908.3275011266899 deg2
#area_obs_S = frac_SGC * area_S # 613.4012361192889 deg2

rds_file = rds_n_file

ra_name_data = 'ALLW_RA'
dec_name_data = 'ALLW_DEC'

ra_name_rds = 'RA'
dec_name_rds = 'DEC'

z_min = float(sys.argv[1])
z_max = float(sys.argv[2])

hduD     = fits.open(data_file)
ra_data    = hduD[1].data[ra_name_data]
dec_data    = hduD[1].data[dec_name_data]
z_data = np.zeros_like(ra_data)
ratelim_data    = hduD[1].data['MASK_2RXS_RATELIM']

#print( hduD[1].data.dtype)

hduR     = fits.open(rds_file)
ra_rds    = hduR[1].data[ra_name_rds]
dec_rds    = hduR[1].data[dec_name_rds]
ratelim_rds    = hduR[1].data['MASK_2RXS_RATELIM']
down_samp = (np.random.random(len(ra_rds))<0.01)

z_s = np.arange(0,2,0.001)
dl_lim = cosmo.luminosity_distance(z_s).to(u.cm)

def write_ascii_clusterings_catalog(out_file, ratelim_min = 0.01, hemisphere='N'):
  spec_SDSS = (hduD[1].data['DR14_Z']>0)&(hduD[1].data['DR14_Z_ERR']>0)&(hduD[1].data['DR14_ZWARNING']==0)
  spec_veron = (hduD[1].data['z_veron']>0)
  z_data[spec_veron] = hduD[1].data['z_veron'][spec_veron]
  z_data[spec_SDSS] = hduD[1].data['DR14_Z'][spec_SDSS]
  spectro = (z_data>z_min)&(z_data<z_max)
  starsD = (hduD[1].data['p_any']>0.5)&(hduD[1].data['2RXS_ExiML']>10)& (hduD[1].data['mask_Tycho20Vmag10']==False)& (hduD[1].data['mask_Tycho210Vmag11']==False)& (hduD[1].data['mask_bright_object_rykoff']==False) & (z_data>0) & (z_data<z_max)
  starsR = (hduR[1].data['mask_Tycho20Vmag10']==False) &(hduR[1].data['mask_Tycho210Vmag11']==False) &(hduR[1].data['mask_bright_object_rykoff']==False)  
 
  if hemisphere=='N':
    sel_data = (z_data<z_max)&(ratelim_data>ratelim_min) &(hduD[1].data['mask_dr14_N'])  &(hduD[1].data['mask_dr14_N_w']>0.9) &(spectro)&(starsD) 
    sel_rds = (down_samp)&(ratelim_rds>ratelim_min) &(hduR[1].data['mask_dr14_N']) &(hduR[1].data['mask_dr14_N_w']>0.9) &(starsR) 

  if hemisphere=='S':
    sel_data = (z_data<z_max)&(ratelim_data>ratelim_min) &(hduD[1].data['mask_dr14_S'])  &(hduD[1].data['mask_dr14_S_w']>0.9) &(spectro)&(starsD)
    sel_rds = (down_samp)&(ratelim_rds>ratelim_min) &(hduR[1].data['mask_dr14_S']) &(hduR[1].data['mask_dr14_S_w']>0.9) &(starsR) 

  # ra figure
  N_data = len(ra_data[sel_data]) 
  N_rds = len(ra_rds[sel_rds]) 
  print("D,R=",N_data, N_rds)
  dz=0.05
  zs=np.arange(0., z_max + dz, dz)
  nn,bb = np.histogram(z_data[sel_data], bins=zs)#, weights=1./w_col.array)
  nz=interp1d((zs[1:]+zs[:-1])/2.,nn)
  rdsz=[]
  for i in range(1,len(zs)-1,1):
    inter=np.random.uniform(low=zs[i]-dz/2., high=zs[i]+dz/2., size=int( 1000* nz( zs[i] )))
    rdsz.append(inter)

  rds=np.hstack((rdsz))
  np.random.shuffle(rds)
  RR=rds[:N_rds]#-dz/2.
  print("RR=",len(rds), len(RR))

  np.savetxt(out_file+str(z_min)+"_"+str(z_max)+"_"+str(ratelim_min)+"_"+hemisphere+".data", np.transpose([ ra_data[sel_data], dec_data[sel_data], z_data[sel_data], np.ones_like(z_data[sel_data]) ]))
  np.savetxt(out_file+str(z_min)+"_"+str(z_max)+"_"+str(ratelim_min)+"_"+hemisphere+".random", np.transpose([ ra_rds[sel_rds], dec_rds[sel_rds], RR, np.ones_like(RR) ]))

print(' / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / ')
print('NORTH')
out_file  = join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/clustering_catalogs/clustering_agn_N_RL_')
write_ascii_clusterings_catalog(out_file, ratelim_min = 0.005, hemisphere='N')
print(' / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / / / / / / / / / // / / / / / / / / / / ')
print('SOUTH')
rds_file = rds_s_file

hduR     = fits.open(rds_file)
ra_rds    = hduR[1].data[ra_name_rds]
dec_rds    = hduR[1].data[dec_name_rds]
ratelim_rds    = hduR[1].data['MASK_2RXS_RATELIM']
down_samp = (np.random.random(len(ra_rds))<0.01)
write_ascii_clusterings_catalog(out_file, ratelim_min = 0.005, hemisphere='S')


