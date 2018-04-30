import numpy as n
import glob

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

import os
import time
import sys
import astropy.io.fits as fits
import astropy.cosmology as co
import astropy.units as u
cosmo = co.Planck13

import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()

from scipy.stats import scoreatpercentile
from scipy.interpolate import interp1d
#out_dir = os.path.join(os.environ['HOME'], "wwwDir", "eRoMok", "AGN", "HG_SMF")

import astropy.io.fits as fits

z_min = 0.# float(sys.argv[1]) # 1.
z_max = 0.3#float(sys.argv[2]) #1.5
#LX_min = 4e42
percentage_allowed = 25. # float(sys.argv[3]) # 20.
tiling_efficiency = 0.95

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

area_obs_N = 908.3275011266899 #deg2
area_obs_S = 613.4012361192889 #deg2
sky_fraction = (area_obs_N+area_obs_S)*n.pi/129600.
print("sky fraction", sky_fraction)
#0.05345070834232634 # data['sky_fraction'][0]

#python plot_Xray_AGN_model.py 0.001 0.2 10.

path_2_cat = os.path.join(os.environ['OBS_REPO'],'SDSS/dr14/spiders/target', "2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask_DR14area_SPM_QSFIT_DR14Q.fits" )

#path_2_cat = os.path.join("/home/comparat/data/spiders", "2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask_DR14area_SPM_QSFIT_DR14Q.fits" )

data = fits.open(path_2_cat)[1].data

DL = cosmo.luminosity_distance(data['DR14_Z']).to(u.cm).value

LX = data['2RXS_SRC_FLUX'] * DL**2 *4. * n.pi
LX_min = 10**(-12.7)* cosmo.luminosity_distance(z_max).to(u.cm).value**2. *4. * n.pi

mass_err = (data['Salpeter_MILES_UVextended_stellar_mass_up'] - data['Salpeter_MILES_UVextended_stellar_mass_low'] )/data['Salpeter_MILES_UVextended_stellar_mass']
mass = data['Salpeter_MILES_UVextended_stellar_mass']

# In a redshift bin
# for a minimum LX
z_bin = (data['DR14_Z']>z_min)&(data['DR14_Z']<z_max)&(LX>LX_min)&(data['mask_dr14_N'])

#agn_vi_class = (z_bin)&(mass>0)&(mass_err>0)&(data['DR14_Z']>0) # ((data['VI_CLASS']=='NLAGN') | (data['VI_CLASS']=='GAL'))&


sel_data_N = (data['mask_dr14_N'])  &(data['mask_dr14_N_w']>0.9) 

sel_data_S = (data['mask_dr14_S'])  &(data['mask_dr14_S_w']>0.9) 

all_sky = (z_bin) & (data['p_any']>0.5) & (data['2RXS_ExiML']>10) & (data['mask_Tycho20Vmag10']==False) & (data['mask_Tycho210Vmag11']==False) & (data['mask_bright_object_rykoff']==False) & ( (sel_data_N) | (sel_data_S) ) 

agn_all = (all_sky)&(mass>0)&(mass_err>0)&(data['DR14_Z']>0)&(mass_err!=n.nan)&(mass!=n.nan)

max_mass_err = scoreatpercentile(mass_err[agn_all], percentage_allowed)

agn = (agn_all) & (all_sky) & (mass_err<max_mass_err)

frac_error_FX = n.median(data['2RXS_SRC_FLUX_ERR'][agn])/n.median(data['2RXS_SRC_FLUX'][agn])
LX_min_up  = LX_min * (1+frac_error_FX)
LX_min_low = LX_min * (1-frac_error_FX)

fraction_with_mass = 1. * len(agn.nonzero()[0]) / len(all_sky.nonzero()[0]) 

print(len(data['DR14_Z'][all_sky]), "targets with Z", len(data['DR14_Z'][agn]), "AGN with mass")
zMean = n.mean(data['DR14_Z'][agn])

volume = (cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min)).value

for pl, mj, fib in zip(data['DR14_PLATE'][agn], data['DR14_MJD'][agn], data['DR14_FIBERID'][agn]):
	os.system("cp ~/SDSS/agn_spiders/*"+str(pl).zfill(4)+"-"+str(mj).zfill(5)+"-"+str(fib).zfill(4)+"* ~/wwwDir/eRoMok/agn_spectra/")