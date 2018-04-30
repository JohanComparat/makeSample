import os
import sys
from os.path import join
import numpy as n
import matplotlib.pyplot as p
from scipy.interpolate import interp1d

out_dir  = join(os.environ['HOME'], 'wwwDir/eRoMok/clustering/data/')
path_2_FX_N = os.path.join( out_dir, 'clustering_agn__N0.0_0.35_0.005_N.2RXS_SRC_FLUX.histogram')
path_2_FX_S = os.path.join( out_dir, 'clustering_agn__S0.0_0.35_0.005_S.2RXS_SRC_FLUX.histogram')
path_2_NZ_N = os.path.join( out_dir, 'clustering_agn__N0.0_0.35_0.005_N.redshift.histogram')
path_2_NZ_S = os.path.join( out_dir, 'clustering_agn__S0.0_0.35_0.005_S.redshift.histogram')


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

import astropy.units as u
import astropy.cosmology as cc
cosmo = cc.Planck13

import lc_setup
import h5py    # HDF5 support
path_2_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_C1.hdf5'
print('opens data')


dz=0.05
z_min=0. 
z_max=0.36
zs=n.arange(z_min, z_max + dz, dz)

f = h5py.File(path_2_lc, 'r')

p.figure(1, (5,5))
p.axes([0.17,0.17,0.78,0.78])

# the mock catalog, FX>-12.52
is_gal = (f['/sky_position/selection'].value) & (f['/agn_properties/agn_activity'].value==1) & (f['/agn_properties/rxay_flux_05_20'].value>10**(-12.52))# &(f['/sky_position/redshift_S'].value>zmin)&(f['/sky_position/redshift_S'].value<zmax)
model_z=f['/sky_position/redshift_S'].value[is_gal]
NN = n.histogram(model_z, bins=zs)[0]
density = NN/3044.
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN**(-0.5), xerr=dz/2., label='mock, FX>-12.52')

# the mock catalog, FX>-127
is_gal = (f['/sky_position/selection'].value) & (f['/agn_properties/agn_activity'].value==1) & (f['/agn_properties/rxay_flux_05_20'].value>10**(-12.7))# &(f['/sky_position/redshift_S'].value>zmin)&(f['/sky_position/redshift_S'].value<zmax)
model_z=f['/sky_position/redshift_S'].value[is_gal]
NN = n.histogram(model_z, bins=zs)[0]
density = NN/3044.
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN**(-0.5), xerr=dz/2., label='mock, FX>-12.7')

z0, z1, density, NN = n.loadtxt(path_2_NZ_N, unpack=True)
p.errorbar((z0+z1)/2., density, yerr=density*NN**(-0.5), xerr=dz/2., label='SPIDERS NGC')
z0, z1, density, NN = n.loadtxt(path_2_NZ_S, unpack=True)
p.errorbar((z0+z1)/2., density, yerr=density*NN**(-0.5), xerr=dz/2., label='SPIDERS SGC')
p.xlabel('redshift')
p.ylabel('N / deg2,  dz=0.05 ')
p.yscale('log')
p.ylim((1e-4, 1))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( out_dir,"histogram_redshift.png"))
p.clf()

dlogF = 0.2
fx_bins = n.arange(-14,-10, dlogF)

p.figure(1, (5,5))
p.axes([0.17,0.17,0.78,0.78])

# the mock catalog, FX>-12.52
is_gal = (f['/sky_position/selection'].value) & (f['/agn_properties/agn_activity'].value==1) & (f['/agn_properties/rxay_flux_05_20'].value>10**(-12.52))# &(f['/sky_position/redshift_S'].value>zmin)&(f['/sky_position/redshift_S'].value<zmax)
flux=f['/agn_properties/rxay_flux_05_20'].value[is_gal]
NN=n.histogram(flux, bins=10**fx_bins)[0]
density = NN/3044.
p.errorbar((fx_bins[1:]+fx_bins[:-1])/2., density, yerr=density*NN**(-0.5), xerr=dlogF/2., label='mock, FX>-12.52')

# the mock catalog, FX>-127
is_gal = (f['/sky_position/selection'].value) & (f['/agn_properties/agn_activity'].value==1) & (f['/agn_properties/rxay_flux_05_20'].value>10**(-12.7))# &(f['/sky_position/redshift_S'].value>zmin)&(f['/sky_position/redshift_S'].value<zmax)
flux=f['/agn_properties/rxay_flux_05_20'].value[is_gal]
NN=n.histogram(flux, bins=10**fx_bins)[0]
density = NN/3044.
p.errorbar((fx_bins[1:]+fx_bins[:-1])/2., density, yerr=density*NN**(-0.5), xerr=dlogF/2., label='mock, FX>-12.7')

z0, z1, density, NN = n.loadtxt(path_2_FX_N, unpack=True)
p.errorbar((z0+z1)/2., density, yerr=density*NN**(-0.5), xerr=dlogF/2., label='SPIDERS NGC')
z0, z1, density, NN = n.loadtxt(path_2_FX_S, unpack=True)
p.errorbar((z0+z1)/2., density, yerr=density*NN**(-0.5), xerr=dlogF/2., label='SPIDERS SGC')
p.xlabel('log(Flux 0.5-2 keV)')
p.ylabel('N / deg2,  dz=0.05 ')
p.yscale('log')
p.grid()
p.ylim((1e-4, 1))
p.xlim((-13.5, -10.))
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( out_dir,"histogram_flux.png"))
p.clf()

