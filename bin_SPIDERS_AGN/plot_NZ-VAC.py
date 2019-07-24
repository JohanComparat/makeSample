import astropy.io.fits as fits
import os, sys, glob
from os.path import join
#import pymangle as mangle
import numpy as np
import matplotlib.pyplot as p
from scipy.interpolate import interp1d

from astropy_healpix import healpy 

import astropy.units as u
import astropy.cosmology as cc
cosmo = cc.Planck13
from astropy.coordinates import SkyCoord

#from sklearnp.neighbors import BallTree, DistanceMetric
from astropy.table import Table,unique,Column
from math import radians, cos, sin, asin, sqrt, pi

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

agn_clustering_dir = '/data36s/comparat/AGN_clustering'
#agn_clustering_dir = '/home/comparat/data/AGN_clustering'

catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )
figure_dir = os.path.join(os.environ['GIT_MAKESAMPLE'], 'figures', 'agn', 'figures_VAC' )
#figure_dir = os.path.join(os.environ['HOME'], 'wwwDir', 'stuff' )

if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)

path_2_rass_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_2RXS_DR16.fits')
hd_rass = fits.open(path_2_rass_cat)[1].data

path_2_xmmsl_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_XMMSL2_DR16.fits')
hd_xmmsl = fits.open(path_2_xmmsl_cat)[1].data


goodZ = (hd_rass['Z_BEST']>-1)
Z_RASS = hd_rass['Z_BEST'][goodZ]

goodZ2 = (hd_xmmsl['Z_BEST']>-1)
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2]


dz=0.1
z_min=0.0005 
z_max=np.max([np.max(Z_RASS),np.max(Z_XMMSL2)])
zs=np.arange(z_min, z_max + dz, dz)
xzs=(zs[1:]+zs[:-1])/2.



def get_z_arrs(path_2_mock):
	hd_mock = fits.open(path_2_mock)[1].data
	fx = hd_mock['agn_FX_soft'] 
	Z_mock_1em12  = np.histogram( hd_mock['redshift_R'][fx>10**(-12.4)] , bins=zs)[0]
	Z_mock_1em125 = np.histogram( hd_mock['redshift_R'][fx>10**(-12.6)] , bins=zs)[0]
	Z_mock_1em127 = np.histogram( hd_mock['redshift_R'][fx>10**(-12.8)] , bins=zs)[0]
	Z_mock_1em13  = np.histogram( hd_mock['redshift_R'][fx>10**(-13.0)] , bins=zs)[0]
	Z_mock_1em14  = np.histogram( hd_mock['redshift_R'][fx>10**(-13.2)] , bins=zs)[0]
	return Z_mock_1em12, Z_mock_1em125, Z_mock_1em127, Z_mock_1em13, Z_mock_1em14

sim_list=np.array(glob.glob(os.path.join(os.environ['MD10'],'cat_AGN_all/00035?.fit')))
simulation_area = healpy.nside2pixarea(8, degrees=True) * len(sim_list)

ZS_data = np.sum(np.array([get_z_arrs(p2s) for p2s in sim_list]), axis=0)/simulation_area

p.figure(0, (10,5))
p.axes([0.12, 0.12, 0.78, 0.78])

p.plot(xzs, ZS_data[0], ls='dashed', lw=2, label='mock, FX>-12.4')
p.plot(xzs, ZS_data[1], ls='dashed', lw=2, label='mock, FX>-12.6')
p.plot(xzs, ZS_data[2], ls='dashed', lw=2, label='mock, FX>-12.8')
p.plot(xzs, ZS_data[3], ls='dotted', lw=2, label='mock, FX>-13.0')
p.plot(xzs, ZS_data[4], ls='dotted', lw=2, label='mock, FX>-13.2')

NN_2rxs = np.histogram(Z_RASS, bins=zs)[0]
density = NN_2rxs/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs**(-0.5), xerr=dz/2., label='2RXS', fmt='none')

NN_xmmxxl = np.histogram(Z_XMMSL2, bins=zs)[0]
density = NN_xmmxxl/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_xmmxxl**(-0.5), xerr=dz/2., label='XMMSL2', fmt='none')

p.xlabel('redshift')
p.ylabel('N/deg2  dz=0.2 ')
p.yscale('log')
p.xlim((0,4))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( figure_dir, "histogram_redshift.png"))
p.clf()

np.savetxt(os.path.join( figure_dir, "histogram_redshift.txt"), np.transpose([zs[:-1], zs[1:], NN_xmmxxl, NN_2rxs]), delimiter=" & ", newline=" \\\\ \n", fmt='%10.1f')
