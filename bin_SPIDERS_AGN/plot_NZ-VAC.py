"""
python 
"""
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
env = sys.argv[1]
catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )
figure_dir = os.path.join(os.environ['GIT_MAKESAMPLE'], 'figures', 'agn', 'figures_VAC' )
#figure_dir = os.path.join(os.environ['HOME'], 'wwwDir', 'stuff' )

if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)

path_2_rass_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_2RXS_DR16.fits')
hd_rass = fits.open(path_2_rass_cat)[1].data

path_2_xmmsl_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_XMMSL2_DR16.fits')
hd_xmmsl = fits.open(path_2_xmmsl_cat)[1].data


def get_selection(data):
	goodZ = (data['Z_BEST']>-1)&(data['CONF_BEST']==3)
	bl  = (goodZ) & ( ( ( data['CLASS_BEST']=='BLAGN') & (data['DR16_CLASS']=='QSO') ) | ( data['CLASS_BEST']=='QSO') )
	nl  = (goodZ) & ( ( data['CLASS_BEST']=='NLAGN') | ( data['CLASS_BEST']=='GALAXY') | (( data['CLASS_BEST']=='BLAGN') & (data['DR16_CLASS']=='GALAXY')) ) 
	ar  = (goodZ) & ( ( data['CLASS_BEST']=='BLAZAR') | ( data['CLASS_BEST']=='BLLAC') )
	bal = (goodZ) & ( ( data['CLASS_BEST']=='BALQSO') | ( data['CLASS_BEST']=='QSO_BAL') )
	no  = (goodZ) & (data['CLASS_BEST']=='NONE')
	st  =(goodZ) & (data['CLASS_BEST']=='STAR')
	print(len(data[goodZ]))
	print(len(data[ bl | nl | ar | bal | no | st ]))
	print('bl ', len(data[bl]))
	print('nl ', len(data[nl]))
	print('ar ', len(data[ar]))
	print('bal', len(data[bal]))
	print('no ', len(data[no]))
	print('st ', len(data[st]))
	return goodZ, bl, nl, ar, bal, no, st

goodZ, bl, nl, ar, bal, no, st = get_selection(hd_rass)

goodZ2, bl2, nl2, ar2, bal2, no2, st2 = get_selection(hd_xmmsl)

Z_RASS = hd_rass['Z_BEST'][goodZ]
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2]

dz=0.1
z_min=0.0005 
z_max=np.max([np.max(Z_RASS),np.max(Z_XMMSL2)])
zs=np.arange(z_min, z_max + dz, dz)
xzs=(zs[1:]+zs[:-1])/2.

def get_z_arrs(path_2_mock):
	print(path_2_mock)
	hd_mock = fits.open(path_2_mock)[1].data
	fx = hd_mock['agn_FX_soft'] 
	mag = hd_mock['AGN_SDSS_r_magnitude']
	AGN_type = hd_mock['AGN_type']
	t1 = (AGN_type==11)|(AGN_type==12) # optically un-obscured
	t2 = (AGN_type==21)|(AGN_type==22) # optically obscured
	Z_mock_1em125_a = np.histogram( hd_mock['redshift_R'][(fx>10**(-12.5)) & (mag>15) & (mag<21.)] , bins=zs)[0]
	Z_mock_1em125_1 = np.histogram( hd_mock['redshift_R'][(fx>10**(-12.5)) & (mag>15) & (mag<21.) & (t1)] , bins=zs)[0]
	Z_mock_1em125_2 = np.histogram( hd_mock['redshift_R'][(fx>10**(-12.5)) & (mag>15) & (mag<21.) & (t2)] , bins=zs)[0]
	Z_mock_1em13_a  = np.histogram( hd_mock['redshift_R'][(fx>10**(-13.0)) & (mag>15) & (mag<21.)] , bins=zs)[0]
	Z_mock_1em13_1  = np.histogram( hd_mock['redshift_R'][(fx>10**(-13.0)) & (mag>15) & (mag<21.) & (t1)] , bins=zs)[0]
	Z_mock_1em13_2  = np.histogram( hd_mock['redshift_R'][(fx>10**(-13.0)) & (mag>15) & (mag<21.) & (t2)] , bins=zs)[0]
	return Z_mock_1em125_a, Z_mock_1em125_1, Z_mock_1em125_2, Z_mock_1em13_a, Z_mock_1em13_1, Z_mock_1em13_2

sim_list5=np.array(glob.glob(os.path.join(os.environ[env],'cat_AGN_all/0003??.fit')))
sim_list6=np.array(glob.glob(os.path.join(os.environ[env],'cat_AGN_all/0004??.fit')))
sim_list = np.hstack((sim_list5,sim_list6))
simulation_area = healpy.nside2pixarea(8, degrees=True) * len(sim_list)

ZS_data = np.sum(np.array([get_z_arrs(p2s) for p2s in sim_list]), axis=0)/simulation_area

p.figure(0, (10,5))
p.axes([0.12, 0.12, 0.78, 0.78])

p.plot(xzs, ZS_data[0], ls='solid', lw=2, label='m, FX>-12.5')
p.plot(xzs, ZS_data[1], ls='dashed', lw=2, label='m, FX>-12.5 t1')
p.plot(xzs, ZS_data[2], ls='dotted', lw=2, label='m, FX>-12.5 t2')
p.plot(xzs, ZS_data[3], ls='solid', lw=2, label='m, FX>-13')
p.plot(xzs, ZS_data[4], ls='dashed', lw=2, label='m, FX>-13 t1')
p.plot(xzs, ZS_data[5], ls='dotted', lw=2, label='m, FX>-13 t2')
#
Z_RASS = hd_rass['Z_BEST'][goodZ]
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2]
#
NN_2rxs = np.histogram(Z_RASS, bins=zs)[0]
density = NN_2rxs/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs**(-0.5), xerr=dz/2., label='2RXS', fmt='none')
#
NN_xmmxxl = np.histogram(Z_XMMSL2, bins=zs)[0]
density = NN_xmmxxl/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_xmmxxl**(-0.5), xerr=dz/2., label='XMMSL2', fmt='none')

# Type 1
Z_RASS = hd_rass['Z_BEST'][goodZ & bl]
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2 & bl2]
#
NN_2rxs_t1 = np.histogram(Z_RASS, bins=zs)[0]
density = NN_2rxs_t1/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs_t1**(-0.5), xerr=dz/2., label='2RXS t1', fmt='none')
#
NN_xmmxxl_t1 = np.histogram(Z_XMMSL2, bins=zs)[0]
density = NN_xmmxxl_t1/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_xmmxxl_t1**(-0.5), xerr=dz/2., label='XMMSL2 t1', fmt='none')

# Type 2
Z_RASS = hd_rass['Z_BEST'][goodZ & nl]
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2 & nl2]
#
NN_2rxs_t2 = np.histogram(Z_RASS, bins=zs)[0]
density = NN_2rxs_t2/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs_t2**(-0.5), xerr=dz/2., label='2RXS t2', fmt='none')
#
NN_xmmxxl_t2 = np.histogram(Z_XMMSL2, bins=zs)[0]
density = NN_xmmxxl_t2/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_xmmxxl_t2**(-0.5), xerr=dz/2., label='XMMSL2 t2', fmt='none')

p.xlabel('redshift')
p.ylabel('N/deg2  dz=0.2 ')
p.yscale('log')
p.xlim((0,4))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( figure_dir, env+"_histogram_redshift.png"))
p.clf()


np.savetxt(os.path.join( figure_dir, env+"_histogram_redshift.txt"), np.transpose([zs[:-1], zs[1:], NN_xmmxxl, NN_2rxs]), delimiter=" & ", newline=" \\\\ \n", fmt='%10.1f')

p.figure(0, (10,5))
p.axes([0.12, 0.12, 0.78, 0.78])

p.plot(xzs, ZS_data[3], ls='solid', lw=2, label='m, FX>-13')
p.plot(xzs, ZS_data[4], ls='dashed', lw=2, label='m, FX>-13 t1')
p.plot(xzs, ZS_data[5], ls='dotted', lw=2, label='m, FX>-13 t2')
#
Z_RASS = hd_rass['Z_BEST'][goodZ]
#
NN_2rxs = np.histogram(Z_RASS, bins=zs)[0]
density = NN_2rxs/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs**(-0.5), xerr=dz/2., label='2RXS', fmt='none')

# Type 1
Z_RASS = hd_rass['Z_BEST'][goodZ & bl]
#
NN_2rxs_t1 = np.histogram(Z_RASS, bins=zs)[0]
density = NN_2rxs_t1/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs_t1**(-0.5), xerr=dz/2., label='2RXS t1', fmt='none')

# Type 2
Z_RASS = hd_rass['Z_BEST'][goodZ & nl]
#
NN_2rxs_t2 = np.histogram(Z_RASS, bins=zs)[0]
density = NN_2rxs_t2/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs_t2**(-0.5), xerr=dz/2., label='2RXS t2', fmt='none')

p.xlabel('redshift')
p.ylabel('N/deg2  dz=0.2 ')
p.yscale('log')
p.xlim((0,4))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( figure_dir, env+"_histogram_redshift_2RXS.png"))
p.clf()


p.figure(0, (10,5))
p.axes([0.12, 0.12, 0.78, 0.78])

p.plot(xzs, ZS_data[0], ls='solid', lw=2, label='m, FX>-12.5')
p.plot(xzs, ZS_data[1], ls='dashed', lw=2, label='m, FX>-12.5 t1')
p.plot(xzs, ZS_data[2], ls='dotted', lw=2, label='m, FX>-12.5 t2')
#
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2]
#
NN_xmmxxl = np.histogram(Z_XMMSL2, bins=zs)[0]
density = NN_xmmxxl/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_xmmxxl**(-0.5), xerr=dz/2., label='XMMSL2', fmt='none')

# Type 1
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2 & bl2]
#
NN_xmmxxl_t1 = np.histogram(Z_XMMSL2, bins=zs)[0]
density = NN_xmmxxl_t1/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_xmmxxl_t1**(-0.5), xerr=dz/2., label='XMMSL2 t1', fmt='none')

# Type 2
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2 & nl2]
#
NN_xmmxxl_t2 = np.histogram(Z_XMMSL2, bins=zs)[0]
density = NN_xmmxxl_t2/5128.#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_xmmxxl_t2**(-0.5), xerr=dz/2., label='XMMSL2 t2', fmt='none')

p.xlabel('redshift')
p.ylabel('N/deg2  dz=0.2 ')
p.yscale('log')
p.xlim((0,4))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( figure_dir, env+"_histogram_redshift_XMMSL2.png"))
p.clf()

