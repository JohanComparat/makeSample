"""
python plot_LX_Z-VAC.py 

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
cosmo = cc.Planck15
from astropy.coordinates import SkyCoord

#from sklearnp.neighbors import BallTree, DistanceMetric
from astropy.table import Table,unique,Column
from math import radians, cos, sin, asin, sqrt, pi

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

#agn_clustering_dir = '/data36s/comparat/AGN_clustering'
agn_clustering_dir = '/home/comparat/data/AGN_clustering'
#env = sys.argv[1]
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

DL_RASS = cosmo.luminosity_distance(Z_RASS).to(u.cm)
DL_XMMSL2 = cosmo.luminosity_distance(Z_XMMSL2).to(u.cm)


FX_RASS = hd_rass['RXS_SRC_FLUX'][goodZ]
FX_XMMSL2 = np.zeros_like(hd_xmmsl['XMMSL2_FLUX_B8'][goodZ2])

ok = (hd_xmmsl['XMMSL2_FLUX_B6'][goodZ2]>0)
FX_XMMSL2[ok] = hd_xmmsl['XMMSL2_FLUX_B6'][goodZ2][ok]*1e-12

ok = (hd_xmmsl['XMMSL2_FLUX_B7'][goodZ2]>0)
FX_XMMSL2[ok] = hd_xmmsl['XMMSL2_FLUX_B7'][goodZ2][ok]*1e-12

ok = (hd_xmmsl['XMMSL2_FLUX_B8'][goodZ2]>0)
FX_XMMSL2[ok] = hd_xmmsl['XMMSL2_FLUX_B8'][goodZ2][ok]*1e-12

LX_RASS = FX_RASS * 4*np.pi * DL_RASS**2.
LX_XMMSL2 = FX_XMMSL2 * 4*np.pi * DL_XMMSL2**2.

zs = np.arange(0,4,0.01)
DL_zs = cosmo.luminosity_distance(zs).to(u.cm)
#LX_1e130 = 10**(-13.0)*4*np.pi*DL_zs**2.
#LX_1e125 = 10**(-12.5)*4*np.pi*DL_zs**2.
#LX_1e120 = 10**(-12.0)*4*np.pi*DL_zs**2.

p.figure(0, (10,5))
p.axes([0.12, 0.12, 0.78, 0.78])

p.plot(Z_RASS, LX_RASS, 'k+', label='2RXS')
p.plot(Z_XMMSL2, LX_XMMSL2, 'rx', label='XMMSL2')

p.plot(zs, 10**(-13.0)*4*np.pi*DL_zs**2., ls='solid', label='-13')
p.plot(zs, 10**(-12.5)*4*np.pi*DL_zs**2., ls='solid', label='-12.5')
p.plot(zs, 10**(-12.0)*4*np.pi*DL_zs**2., ls='solid', label='-12')


p.xlabel('redshift')
p.ylabel(r'$\log_{10}(L_X/[erg/s])$')
p.xscale('log')
p.yscale('log')
p.xlim((0.01,4.2))
p.ylim((1e40,1e48))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( figure_dir, "LX_redshift.png"))
p.clf()



