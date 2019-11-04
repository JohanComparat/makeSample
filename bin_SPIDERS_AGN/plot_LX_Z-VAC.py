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


area_eboss_dr16 = 5128.

MIN_SDSS_FIBER2MAG_i = 17.0
MAX_SDSS_FIBER2MAG_i = 22.5
MIN_SDSS_MODELMAG_i  = 16.0

sig_fx = 0.3 #1e-13
sig_fx = 0.2 #1e-13
sig_fx = 0.1 #sig_fx = 2e-13

MAG_MIN_LIM = 15
MAG_LIM = 22.5

#agn_clustering_dir = '/data36s/comparat/AGN_clustering'
agn_clustering_dir = '/home/comparat/data/AGN_clustering'
#env = sys.argv[1]
catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )
figure_dir = os.path.join(os.environ['GIT_MAKESAMPLE'], 'figures', 'agn', 'figures_VAC' )
#figure_dir = os.path.join(os.environ['HOME'], 'wwwDir', 'stuff' )

if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)

path_2_rass_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_2RXS_DR16.fits')
data_2RXS = fits.open(path_2_rass_cat)[1].data

path_2_xmmsl_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_XMMSL2_DR16.fits')
data_XMM = fits.open(path_2_xmmsl_cat)[1].data


def get_arrays_xmmsl2(path_2_cat, EXI_ML_min = 10):
	data_2RXS_i = fits.open(path_2_cat)[1].data
	targets = (data_2RXS_i['NWAY_match_flag']!=2) & (data_2RXS_i['XMMSL2_IN_BOSS']==1) & (data_2RXS_i['FLAG_SDSSv5b_best']==1) & (data_2RXS_i['NWAY_p_any']>=0.01) & (data_2RXS_i['NUM_SDSS']>=1)
	#  NWAY_match_flag!=2 &&RXS_IN_BOSS==1 &&FLAG_SDSSv5b_best==1 &&
	data_2RXS = fits.open(path_2_cat)[1].data[targets]
	#print(np.unique(data_2RXS['CLASS_BEST']))
	data = data_2RXS
	ra = data['XMMSL2_RA']
	dec = data['XMMSL2_DEC']
	ExiML = np.transpose([
		data['XMMSL2_DET_ML_B6'] , 
		data['XMMSL2_DET_ML_B7'] , 
		data['XMMSL2_DET_ML_B8'] ])
	ExiML[np.isnan(ExiML)]=0
	ExiML_arr = np.max(ExiML,axis=1) 
	high_conf = (ExiML_arr >= EXI_ML_min )
	targeted  = (high_conf) & (data['SDSS_FIBER2MAG_i']>=MIN_SDSS_FIBER2MAG_i) & (data['SDSS_FIBER2MAG_i']<=MAX_SDSS_FIBER2MAG_i) & (data['SDSS_MODELMAG_i']>=MIN_SDSS_MODELMAG_i)    
	observed  = (targeted) & (data['DR16_MEMBER']==1)
	c1 = (observed) & ((data['Z_BEST']>-0.5) | ((data['DR16_Z']>-0.5) & (data['DR16_Z_ERR']>0)))
	c2 = (c1) & (data['CONF_BEST']==3)
	c3 = (c1) & (data['CONF_BEST']==2) & ((data['CLASS_BEST']=='BLAZAR')|(data['CLASS_BEST']=='BLLAC'))
	c4 = (c1) & (data['DR16_SN_MEDIAN_ALL']>=2) & (data['DR16_ZWARNING']==0 )
	c5 = (c1) & (data['CONF_BEST']==2) & (data['DR16_ZWARNING']==0 ) & (data['VI_REINSPECT_FLAG'] == 0) & (data['VI_NINSPECTORS']>2)
	c6 = (c1) & (data['CONF_BEST']==2) & (data['DR16_ZWARNING']==0 ) & (data['VI_AM_RECONCILED']==1)
	idZ =  (c2) | (c3) | (c4) | (c5) | (c6) 
	blazars_noZ = (idZ) & ((data['CLASS_BEST']=='BLAZAR')|(data['CLASS_BEST']=='BLLAC')) & (data['CONF_BEST']<3)
	goodZ = (idZ) & (blazars_noZ == False)
	clusters = (goodZ) & (data['REDMAPPER_Separation']<60) & (abs(data['DR16_Z'] - data['REDMAPPER_Z_LAMBDA'])<0.01)
	stars = (goodZ) & (data['CLASS_BEST']=='STAR')
	agnZ = (goodZ) & (clusters==False) &  (stars==False)
	return data, ra, dec, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, data_2RXS_i


def get_arrays(data = data_2RXS, EXI_ML_min = 6.5):
	#ra = data['RXS_RAJ2000']
	#dec = data['RXS_DEJ2000']
	high_conf = (data['RXS_ExiML'] >= EXI_ML_min )
	targeted  = (high_conf) & (data['SDSS_FIBER2MAG_i']>=MIN_SDSS_FIBER2MAG_i) & (data['SDSS_FIBER2MAG_i']<=MAX_SDSS_FIBER2MAG_i) & (data['SDSS_MODELMAG_i']>=MIN_SDSS_MODELMAG_i)    
	observed  = (targeted) & (data['DR16_MEMBER']==1)
	c1 = (observed) & ((data['Z_BEST']>-0.5) | ((data['DR16_Z']>-0.5) & (data['DR16_Z_ERR']>0)))
	c2 = (c1) & (data['CONF_BEST']==3)
	c3 = (c1) & (data['CONF_BEST']==2) & ((data['CLASS_BEST']=='BLAZAR')|(data['CLASS_BEST']=='BLLAC'))
	c4 = (c1) & (data['DR16_SN_MEDIAN_ALL']>=2) & (data['DR16_ZWARNING']==0 )
	c5 = (c1) & (data['CONF_BEST']==2) & (data['DR16_ZWARNING']==0 ) & (data['VI_REINSPECT_FLAG'] == 0) & (data['VI_NINSPECTORS']>2)
	c6 = (c1) & (data['CONF_BEST']==2) & (data['DR16_ZWARNING']==0 ) & (data['VI_AM_RECONCILED']==1)
	idZ =  (c2) | (c3) | (c4) | (c5) | (c6) 
	blazars_noZ = (idZ) & ((data['CLASS_BEST']=='BLAZAR')|(data['CLASS_BEST']=='BLLAC')) & (data['CONF_BEST']<3)
	goodZ = (idZ) & (blazars_noZ == False)
	clusters = (goodZ) & (data['REDMAPPER_Separation']<60) & (abs(data['DR16_Z'] - data['REDMAPPER_Z_LAMBDA'])<0.01)
	stars = (data['CLASS_BEST']=='STAR')
	agnZ = (goodZ) & (clusters==False) &  (stars==False)
	bl  = (agnZ) & ( ( ( data['CLASS_BEST']=='BLAGN') & (data['DR16_CLASS']=='QSO') ) | ( data['CLASS_BEST']=='QSO') )
	nl  = (agnZ) & ( bl==False)
	#print(len(ra), len(ra[high_conf]) )
	return data, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, bl, nl


path_2_cat = os.path.join(catalog_dir, '2RXS', 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER.fits')
data_2RXS_i = fits.open(path_2_cat)[1].data
targets = (data_2RXS_i['NWAY_match_flag']!=2) & (data_2RXS_i['RXS_IN_BOSS']==1) & (data_2RXS_i['FLAG_SDSSv5b_best']==1) & (data_2RXS_i['NWAY_p_any']>=0.01) & (data_2RXS_i['NUM_SDSS']>=1)
data_2RXS = data_2RXS_i[targets]

data_2RXS, targeted_2RXS, observed_2RXS, goodZ_2RXS, idZ_2RXS, agnZ_2RXS, clusterZ_2RXS, blazarNOZ_2RXS, stars_2RXS, bl_2RXS, nl_2RXS = get_arrays(data_2RXS)


path_2_cat_xmmsl = os.path.join(catalog_dir, 'XMMSL2',  'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_REDMAPPER.fits')

data_XMM, ra_XMM, dec_XMM, targeted_XMM, observed_XMM, goodZ_XMM, idZ_XMM, agnZ_XMM, clusterZ_XMM, blazarNOZ_XMM, stars_XMM, all_data_XMM = get_arrays_xmmsl2(path_2_cat_xmmsl)

Z_RASS = data_2RXS['Z_BEST'][agnZ_2RXS]
Z_XMMSL2 = data_XMM['Z_BEST'][agnZ_XMM]

DL_RASS = cosmo.luminosity_distance(Z_RASS).to(u.cm)
DL_XMMSL2 = cosmo.luminosity_distance(Z_XMMSL2).to(u.cm)

FX_RASS = data_2RXS['RXS_SRC_FLUX'][agnZ_2RXS]
FX_XMMSL2 = np.zeros_like(data_XMM['XMMSL2_FLUX_B8'][agnZ_XMM])

ok = (data_XMM['XMMSL2_FLUX_B6'][agnZ_XMM]>0)
FX_XMMSL2[ok] = data_XMM['XMMSL2_FLUX_B6'][agnZ_XMM][ok]*1e-12

ok = (data_XMM['XMMSL2_FLUX_B7'][agnZ_XMM]>0)
FX_XMMSL2[ok] = data_XMM['XMMSL2_FLUX_B7'][agnZ_XMM][ok]*1e-12

ok = (data_XMM['XMMSL2_FLUX_B8'][agnZ_XMM]>0)
FX_XMMSL2[ok] = data_XMM['XMMSL2_FLUX_B8'][agnZ_XMM][ok]*1e-12

LX_RASS = FX_RASS * 4*np.pi * DL_RASS**2.
LX_XMMSL2 = FX_XMMSL2 * 4*np.pi * DL_XMMSL2**2.

zs = np.arange(0,4,0.01)
DL_zs = cosmo.luminosity_distance(zs).to(u.cm)
#LX_1e130 = 10**(-13.0)*4*np.pi*DL_zs**2.
#LX_1e125 = 10**(-12.5)*4*np.pi*DL_zs**2.
#LX_1e120 = 10**(-12.0)*4*np.pi*DL_zs**2.

p.figure(0, (5.5,5.5))
p.axes([0.15, 0.12, 0.78, 0.78])

#p.plot(Z_RASS, LX_RASS, 'k+', label='2RXS')
p.plot(Z_XMMSL2, LX_XMMSL2, 'kx', label='XMMSL2')

p.plot(zs, 10**(-13.0)*4*np.pi*DL_zs**2., ls='solid', label='-13')
p.plot(zs, 10**(-12.5)*4*np.pi*DL_zs**2., ls='solid', label='-12.5')
p.plot(zs, 10**(-12.0)*4*np.pi*DL_zs**2., ls='solid', label='-12')

p.xlabel('redshift')
p.ylabel(r'$\log_{10}(L_X/[erg/s])$')
p.xscale('log')
p.yscale('log')
xtk = np.array([0.03, 0.06, 0.1, 0.3, 0.06, 1, 2, 4])
p.xticks(xtk, xtk)
p.xlim((0.02,4.2))
p.ylim((1e41,1e47))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( figure_dir, "LX_redshift_XMMSL2.png"))
p.clf()


p.figure(0, (5.5,5.5))
p.axes([0.15, 0.12, 0.78, 0.78])

p.plot(Z_RASS, LX_RASS, 'k+', label='2RXS')
#p.plot(Z_XMMSL2, LX_XMMSL2, 'rx', label='XMMSL2')

p.plot(zs, 10**(-13.0)*4*np.pi*DL_zs**2., ls='solid', label='-13')
p.plot(zs, 10**(-12.5)*4*np.pi*DL_zs**2., ls='solid', label='-12.5')
p.plot(zs, 10**(-12.0)*4*np.pi*DL_zs**2., ls='solid', label='-12')

p.xlabel('redshift')
p.ylabel(r'$\log_{10}(L_X/[erg/s])$')
p.xscale('log')
p.yscale('log')
p.xlim((0.02,4.2))
xtk = np.array([0.03, 0.06, 0.1, 0.3, 0.06, 1, 2, 4])
p.xticks(xtk, xtk)
p.ylim((1e41,1e47))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( figure_dir, "LX_redshift_2RXS.png"))
p.clf()



