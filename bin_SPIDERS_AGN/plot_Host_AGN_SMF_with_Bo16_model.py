import numpy as n
import glob

import h5py    # HDF5 support

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

import os
import time
import sys
import astropy.io.fits as fits
import astropy.cosmology as co
import astropy.units as u
cosmo = co.Planck15

import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()

from scipy.stats import scoreatpercentile
from scipy.interpolate import interp1d
out_dir = os.path.join(os.environ['HOME'], "wwwDir", "eRoMok", "AGN", "HG_SMF")

import astropy.io.fits as fits

z_min = float(sys.argv[1]) # 1.
z_max = float(sys.argv[2]) #1.5
#LX_min = 4e42
percentage_allowed = float(sys.argv[3]) # 20.
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

path_2_cat = os.path.join("/data36s/comparat/SDSS/dr14/spiders/target", "2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask_DR14area_SPM.fits" )

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

nc,bc,pc = p.hist(n.log10(mass[agn]), bins = 1000, cumulative=True)
p.clf()
it_4_bins = interp1d(nc,(bc[1:]+bc[:-1])/2.)

N_per_bin = 16
mbins = n.hstack(( bc[0], it_4_bins(n.arange(1, len(data['DR14_Z'][agn]), N_per_bin)), bc[-1] ))
dlog10M = mbins[1:]-mbins[:-1]
xb = (mbins[1:] + mbins[:-1]) / 2.

NN, bb = n.histogram(n.log10(mass[agn]), bins = mbins)

y = NN / (sky_fraction * volume * fraction_with_mass * tiling_efficiency) / dlog10M
y_err = NN**(-0.5) * y 


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


all_z = n.arange(0, 6, 0.1)

#p.figure(1, (6,6))
#p.axes([0.17, 0.17, 0.78, 0.78])
#p.plot(all_z, xr.f_z(all_z))
#p.xlabel('redshift')
#p.ylabel(r'$f(z)$')
##p.xlim((9., 12.2))
##p.ylim((-9,-2))
#p.grid()
##p.legend(loc=0, frameon=False)
#p.savefig(os.path.join(out_dir, "B016_fz.png"))
#p.clf()

many_MS = n.arange(9.,12.5,0.05)

##################################################3
##################################################3
##################################################3
# SIMUATION COUNTERPART
##################################################3
##################################################3

f = h5py.File('/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3_.hdf5', 'r')

is_gal = (f['/sky_position/selection'].value)
is_agn = (f['/sky_position/selection'].value)&(f['/agn_properties/agn_activity'].value==1)

n_gal = len(f['/sky_position/redshift_S'].value[is_gal])

n_agn = len(f['/sky_position/redshift_S'].value[is_agn])

sel = (is_agn)&(f['/sky_position/redshift_S'].value>z_min)&(f['/sky_position/redshift_S'].value<z_max)&(n.log10(f['/moster_2013_data/stellar_mass'].value)+f['/agn_properties/log_lambda_sar'].value>lxmin)

n.savetxt(out_name, n.transpose([f['/sky_position/RA'].value[sel], f['/sky_position/DEC'].value[sel], f['/sky_position/redshift_S'].value[sel], n.ones_like(f['/sky_position/redshift_S'].value[sel])]) )
print(zmax, lxmin, len(f['/sky_position/RA'].value[sel]))

area = 6.7529257176359*2. * 2* 8.269819492449505
volume_sim = volume * area *n.pi / 129600. 

logM = n.log10(f['/moster_2013_data/stellar_mass'].value[sel])


NN_sim = n.histogram(logM, bins = mbins)[0]
y_sim = NN_sim / volume_sim / dlog10M
y_sim_err = NN_sim**(-0.5) * y_sim 


p.figure(1, (6,6))
p.axes([0.17, 0.17, 0.78, 0.78])
p.plot(data['DR14_Z'], LX, 'ko', label='all')
#p.plot(data['DR14_Z'][agn_all], LX[agn_all], 'r+', label='nlagn or gal')
p.plot(data['DR14_Z'][agn], LX[agn], 'r+')#, label=str(max_mass_err))
p.plot(n.arange(0,2,0.005), 10**(-12.7)* cosmo.luminosity_distance(n.arange(0,2,0.005)).to(u.cm).value**2. *4. * n.pi, label='log Flux>-12.7')

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


overall = n.array([xr.Phi_stellar_mass(MS, zMean) for MS in many_MS])
p.plot(many_MS, overall, label='all')

#ttt = n.array([xr.Phi_stellar_mass_mmin(MS, zMean, llx-MS) for MS in many_MS])
lxlim = n.log10(LX_min_up ) 
ttt_up = n.min(n.array([n.array([xr.Phi_stellar_mass_mmin(MS, zMean, lxlim-MS) for MS in many_MS]), overall]), axis=0) 
lxlim = n.log10(LX_min_low )
ttt_low = n.min(n.array([n.array([xr.Phi_stellar_mass_mmin(MS, zMean, lxlim-MS) for MS in many_MS]), overall]), axis=0) 
lxlim = n.log10(LX_min)
print("LX limit=",lxlim)

p.fill_between(many_MS, y1=ttt_low, y2=ttt_up, label='LX>'+str(n.round(lxlim,2)), alpha=0.5, color='magenta')

p.fill_between(xb, y1=y_sim-y_sim_err, y2=y_sim+y_sim_err, label='mock', alpha=0.5, color='green', alpha=0.5)

p.errorbar(xb, y, xerr=dlog10M/2., yerr=y_err, label=r'spiders $\sigma_M<$'+str(n.round(max_mass_err,3)), fmt=None, color='k')
#p.plot(bb[:-1]+0.1, y*100./percentage_allowed, 'k+', label='spiders 100%')
i_ilb13=0
p.plot(many_MS, smf_ilbert_fun[i_ilb13](10**many_MS), label='Ilb13 ' + str(smf_ilbert_zmin[i_ilb13]) + '<z<' + str(smf_ilbert_zmax[i_ilb13]), ls='dashed' )

def plot_llx(llx, style='dotted'):
  ttt = n.array([xr.Phi_stellar_mass_mmin(MS, zMean, llx-MS) for MS in many_MS])
  new = n.min(n.array([ttt, overall]), axis=0) 
  p.plot(many_MS, new, ls=style, label='LX>'+str(llx))

plot_llx(42)
plot_llx(43)
plot_llx(44)


p.xlabel(r'$log_{10}(M*/M_\odot)$')
p.ylabel(r'$\Phi_*(M) [Mpc^{-3} dlog_{10}M^{-1}]$')
p.xlim((9.5, 12.))
p.ylim((1e-9,1e-2))
p.yscale('log')
p.grid()
p.title('z='+str(n.round(zMean,3)) + ", fraction used=" + str(percentage_allowed) + '%')
p.legend(loc=1, frameon=False)
p.savefig(os.path.join(out_dir, "B016_psi_star_LX_"+str(z_min)+"_"+str(z_max)+"_pc_"+str(percentage_allowed)+".png"))
p.clf()

