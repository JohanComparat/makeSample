import numpy as n
import glob
import matplotlib.pyplot as p
#import h5py
import os
import time
import sys
import astropy.io.fits as fits
import astropy.cosmology as co
cosmo = co.Planck15
import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()
import astropy.units as u
from scipy.stats import scoreatpercentile

#out_dir = os.path.join("/afs/mpe/www/people/comparat/", "eRoMok", "h5", "xray_agn_model" )
out_dir = os.path.join(os.environ['DATA_DIR'],"spiders","agn" )

import astropy.io.fits as fits

z_min = float(sys.argv[1]) # 1.
z_max = float(sys.argv[2]) #1.5
#LX_min = 4e42
percentage_allowed = float(sys.argv[3]) # 20.

fraction_with_mass = 4602./9073.
tiling_efficiency = 0.95

#python plot_Xray_AGN_model.py 0. 0.2 10.

path_2_cat = '/home/comparat/data/spiders/VAC_2RXS_DR14_ALLW_SDSS_GAIA_FIRST_final_rest-frame_with_MS.FITS'

path_2_cat = '/home/comparat/data/spiders/VAC_2RXS_DR14_ALLW_SDSS_GAIA_FIRST_final_SPM.FITS'

data = fits.open(path_2_cat)[1].data

DL = cosmo.luminosity_distance(data['VI_Z']).to(u.cm).value
LX = data['RXS_SRC_FLUX'] * DL**2 *4. * n.pi
LX_min = 2e-13* cosmo.luminosity_distance(z_max).to(u.cm).value**2. *4. * n.pi

mass_err = (data['Chabrier_MILES_stellar_mass_up'] - data['Chabrier_MILES_stellar_mass_low'] )/data['Chabrier_MILES_stellar_mass']
mass = data['Chabrier_MILES_stellar_mass']

# In a redshift bin
# for a minimum LX
z_bin = (data['VI_Z']>z_min)&(data['VI_Z']<z_max)&(LX>LX_min)

agn_vi_class = (z_bin)&((data['VI_CLASS']=='NLAGN') | (data['VI_CLASS']=='GAL'))&(mass>0)&(mass_err>0)&(data['VI_Z']>0)

agn_all = (z_bin)&(mass>0)&(mass_err>0)&(data['VI_Z']>0)

max_mass_err = scoreatpercentile(mass_err[agn_all], percentage_allowed)

agn = (agn_all) & (z_bin) & (mass_err<max_mass_err)

print(len(data['VI_Z'][agn]), "points")
zMean = n.mean(data['VI_Z'][agn])

sky_fraction = 0.05345070834232634 # data['sky_fraction'][0]
volume = (cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min)).value

NN, bb = n.histogram(n.log10(mass[agn]), bins = n.arange(7,12.5,0.25))

y = NN / (sky_fraction * volume * fraction_with_mass * tiling_efficiency)
y_err = NN**(-0.5) * y 


bins = n.arange(6,13,0.1)
xb = (bins[1:] + bins[:-1]) / 2.

smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
ll_dir = os.path.join(os.environ['GIT_NBODY_NPT'], 'data', 'stellar_mass_function')
path_ilbert13_SMF = os.path.join(ll_dir, "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(path_ilbert13_SMF, unpack=True)

smf_ilbert_fun = n.array([
lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
, lambda mass : smf_ilbert13( mass , 10**M_star[1], phi_1s[1]*10**(-3), alpha_1s[1], phi_2s[1]*10**(-3), alpha_2s[1] )
, lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )
, lambda mass : smf_ilbert13( mass , 10**M_star[3], phi_1s[3]*10**(-3), alpha_1s[3], phi_2s[3]*10**(-3), alpha_2s[3] )
, lambda mass : smf_ilbert13( mass , 10**M_star[4], phi_1s[4]*10**(-3), alpha_1s[4], phi_2s[4]*10**(-3), alpha_2s[4] )
, lambda mass : smf_ilbert13( mass , 10**M_star[5], phi_1s[5]*10**(-3), alpha_1s[5], phi_2s[5]*10**(-3), alpha_2s[5] )
, lambda mass : smf_ilbert13( mass , 10**M_star[6], phi_1s[6]*10**(-3), alpha_1s[6], phi_2s[6]*10**(-3), alpha_2s[6] )
, lambda mass : smf_ilbert13( mass , 10**M_star[7], phi_1s[7]*10**(-3), alpha_1s[7], phi_2s[7]*10**(-3), alpha_2s[7] )
])

mbins = n.arange(8,12.5,0.25)

smf_ilbert_zmin = n.array([ 
0.2
, 0.5
, 0.8
, 1.1
, 1.5
, 2.0
, 2.5
, 3.0 ])

smf_ilbert_zmax = n.array([ 
0.5
, 0.8
, 1.1
, 1.5
, 2.0
, 2.5
, 3.0
, 4.0 ])

smf_ilbert_name = n.array([ "Il13 "+str(zmin)+"<z<"+str(zmax) for zmin, zmax in zip(smf_ilbert_zmin,smf_ilbert_zmax) ])


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

all_z = n.arange(0, 6, 0.1)

p.figure(1, (6,6))
p.axes([0.17, 0.17, 0.78, 0.78])
p.plot(all_z, xr.f_z(all_z))
p.xlabel('redshift')
p.ylabel(r'$f(z)$')
#p.xlim((9., 12.2))
#p.ylim((-9,-2))
p.grid()
#p.legend(loc=0, frameon=False)
p.savefig(os.path.join(out_dir, "B016_fz.png"))
p.clf()

many_MS = n.arange(8.7,12.5,0.1)


p.figure(1, (6,6))
p.axes([0.17, 0.17, 0.78, 0.78])
p.plot(data['VI_Z'], LX, 'ko', label='all')
#p.plot(data['VI_Z'][agn_all], LX[agn_all], 'r+', label='nlagn or gal')
p.plot(data['VI_Z'][agn], LX[agn], 'r+')#, label=str(max_mass_err))
p.plot(n.arange(0,2,0.005), 2e-13* cosmo.luminosity_distance(n.arange(0,2,0.005)).to(u.cm).value**2. *4. * n.pi, label='Flux=2e-13')

p.xlabel(r'redshift')
p.ylabel(r'$LX$')
p.xlim((z_min-0.1,z_max+0.1))
p.ylim((1e40,1e47))
p.yscale('log')
p.grid()
p.legend(loc=0, frameon=False)
p.savefig(os.path.join(out_dir, "SPIDERS_AGN_redshift_LX_"+str(z_min)+"_"+str(z_max)+".png"))
p.clf()



p.figure(1, (6,6))
p.axes([0.17, 0.17, 0.78, 0.78])
p.errorbar(bb[:-1]+0.125, y, xerr=0.125, yerr=y_err/2., label=r'spiders $\sigma_M<$'+str(n.round(max_mass_err,3)))
p.plot(bb[:-1]+0.125, y*100./percentage_allowed, 'k+', label='spiders 100%')

overall = n.array([xr.Phi_stellar_mass(MS, zMean) for MS in many_MS])
p.plot(many_MS, overall, label='all')

lxlim = n.log10(LX_min)
new = n.min(n.array([n.array([xr.Phi_stellar_mass_mmin(MS, zMean, lxlim-MS) for MS in many_MS]), overall]), axis=0) 
p.plot(many_MS, new, ls='solid', label='>'+str(n.round(lxlim,2)))

for lxlim in n.arange(42,44.6,1.):
	print(lxlim)
	new = n.min(n.array([n.array([xr.Phi_stellar_mass_mmin(MS, zMean, lxlim-MS) for MS in many_MS]), overall]), axis=0) 
	p.plot(many_MS, new, ls='dashed', label='LX>'+str(lxlim))



p.xlabel(r'$log(M*/M_\odot)$')
p.ylabel(r'$\Phi_*(M)$')
p.xlim((9., 12.5))
p.ylim((1e-9,1e-3))
p.yscale('log')
p.grid()
p.title('z='+str(n.round(zMean,3))+" % most accurate="+str(percentage_allowed))
p.legend(loc=1, frameon=False)
p.savefig(os.path.join(out_dir, "B016_psi_star_LX_"+str(z_min)+"_"+str(z_max)+"_pc_"+str(percentage_allowed)+".png"))
p.clf()

