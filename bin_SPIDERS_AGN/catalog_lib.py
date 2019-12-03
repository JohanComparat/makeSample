"""

library of function to be loaded to retrieve samples from 2RXS and XMMSL2

"""
write_table = False 
update_table = False 

import astropy.io.fits as fits
import os, sys, glob
from os.path import join

import numpy as np

from scipy.interpolate import interp1d
from scipy.stats import norm  

import pymangle 
import time
t0 = time.time()

from astropy.coordinates import SkyCoord

from astropy_healpix import healpy 
import healpy as hp

import astropy.units as u
import astropy.cosmology as cc

from astropy.cosmology import FlatLambdaCDM
#cosmo = cc.Planck13
cosmoMD = FlatLambdaCDM(H0=67.77 * u.km / u.s / u.Mpc,Om0=0.307115)  # , Ob0=0.048206)
h = 0.6777
L_box = 1000.0 / h
cosmo = cosmoMD

print("interpolates z d_comoving and D_L")
z_array = np.arange(0, 7.5, 0.001)
dc_to_z = interp1d(cosmo.comoving_distance(z_array), z_array)
d_L = cosmo.luminosity_distance(z_array)
dl_cm = (d_L.to(u.cm)).value
dL_interpolation = interp1d(z_array, dl_cm)


#from sklearnp.neighbors import BallTree, DistanceMetric
from astropy.table import Table,unique,Column
from math import radians, cos, sin, asin, sqrt, pi


#agn_clustering_dir = '/data36s/comparat/AGN_clustering'
git_dir = os.environ['GIT_MAKESAMPLE']
figure_dir = os.path.join(git_dir, 'figures', 'agn', 'figures_VAC' )
if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)

agn_clustering_dir = '/home/comparat/data/AGN_clustering'

target_dir  = os.path.join(agn_clustering_dir, 'targets'  )
footprint_dir  = os.path.join(agn_clustering_dir, 'footprint'  )
catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )


#path_2_cat_2RXS = os.path.join(catalog_dir, '2RXS', 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER_SPIDERSCODEX_XID.fits')
path_2_cat_2RXS = os.path.join(catalog_dir, '2RXS', 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1.fits')

#path_2_cat_xmmsl = os.path.join(catalog_dir, 'XMMSL2',  'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_REDMAPPER_SPIDERSCODEX_XID.fits')
path_2_cat_xmmsl = os.path.join(catalog_dir, 'XMMSL2',  'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1.fits')

# selection cuts : 
MIN_SDSS_FIBER2MAG_i = 17.0
MAX_SDSS_FIBER2MAG_i = 22.5
MIN_SDSS_MODELMAG_i  = 16.0
EXI_ML_min=6.5

def get_arrays_xmmsl2(path_2_cat_xmmsl = path_2_cat_xmmsl, EXI_ML_min = 6.5, option='none'):
	"""
	Samples and selections in XMMSL2
	"""
	print("============================")
	print("XMMSL2")
	data_XMMSL2_i = fits.open(path_2_cat_xmmsl)[1].data
	targets = (data_XMMSL2_i['NWAY_match_flag']!=2) & (data_XMMSL2_i['XMMSL2_IN_BOSS']==1) & (data_XMMSL2_i['FLAG_SDSSv5b_best']==1) & (data_XMMSL2_i['NWAY_p_any']>=0.01) & (data_XMMSL2_i['NUM_SDSS']>=1)
	#  NWAY_match_flag!=2 &&RXS_IN_BOSS==1 &&FLAG_SDSSv5b_best==1 &&
	data = data_XMMSL2_i[targets]
	#revised = (data['XID_PLATE_BEST']>0)
	#print(data['XID_CLASS_REVISED'][revised])
	#class_best = np.copy( data['CLASS_BEST'] )#.astype('str')
	#class_best[revised] = data['XID_CLASS_REVISED'][revised]#.astype('str')
	#class_best = np.array([ el.strip() for el in class_best]) 
	if option=='full_boss':
		class_best = data['DR16_CLASS'] 
	else:
		class_best = data['merged_class']
	print(np.unique(class_best, return_counts=True))
	#print(np.unique(data_2RXS['CLASS_BEST']))
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
	c3 = (c1) & (data['CONF_BEST']==2) & ((class_best=='BLAZAR')|(class_best=='BLLAC'))
	c4 = (c1) & (data['DR16_SN_MEDIAN_ALL']>=2) & (data['DR16_ZWARNING']==0 )
	c5 = (c1) & (data['CONF_BEST']==2) & (data['DR16_ZWARNING']==0 ) & (data['VI_REINSPECT_FLAG'] == 0) & (data['VI_NINSPECTORS']>2)
	c6 = (c1) & (data['CONF_BEST']==2) & (data['DR16_ZWARNING']==0 ) & (data['VI_AM_RECONCILED']==1)
	idZ =  (c2) | (c3) | (c4) | (c5) | (c6) 
	blazars_noZ = (idZ) & ((class_best=='BLAZAR')|(class_best=='BLLAC')) & (data['CONF_BEST']<3)
	goodZ = (idZ) & (blazars_noZ == False)
	redmapper_cluster = (goodZ) & (data['REDMAPPER_Separation']<60) & (abs(data['DR16_Z'] - data['REDMAPPER_Z_LAMBDA'])<0.01)
	spiders_cluster = (goodZ) & (abs(data['SPIDERSCODEX_SCREEN_CLUZSPEC']-data['Z_BEST'])<0.01)
	clusters = (redmapper_cluster) | (spiders_cluster)
	stars = (goodZ) & (class_best=='STAR')
	agnZ = (goodZ) & (clusters==False) &  (stars==False)
	bl  = (agnZ) & ( (class_best=='BLAGN')  | (class_best=='QSO') )
	other = (agnZ) & ( (class_best=="BLLAC") | ( class_best== "BALQSO")  | ( class_best== "BLAZAR")  | ( class_best== "NONE")  | ( class_best== "QSO_BAL") )  
	nl  = (agnZ) & ( bl==False) & (other==False)
	print(len(data), len(targeted), len(observed), len(goodZ), len(idZ), len(agnZ), len(clusters), len(blazars_noZ), len(stars), len(bl), len(nl) )
	return data, high_conf, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, bl, nl, class_best


def get_arrays_2rxs(path_2_cat_2RXS = path_2_cat_2RXS, EXI_ML_min = 6.5, option='none'):
	"""
	Samples and selections in 2RXS
	"""
	print("============================")
	print("2RXS")
	hd_rass_i = fits.open(path_2_cat_2RXS)[1].data
	targets = (hd_rass_i['NWAY_match_flag']!=2) & (hd_rass_i['RXS_IN_BOSS']==1) & (hd_rass_i['FLAG_SDSSv5b_best']==1) & (hd_rass_i['NWAY_p_any']>=0.01) & (hd_rass_i['NUM_SDSS']>=1)
	data = hd_rass_i[targets]
	# revision of class best by A. Schwope
	#revised = (data['XID_Flux_synt']>0)
	#print(data['XID_CLASS_REVISED'][revised])
	#class_best = np.copy( data['CLASS_BEST'] ).astype('str')
	#class_best[revised] = data['XID_CLASS_REVISED'][revised].astype('str')
	#class_best = np.array([ el.strip() for el in class_best]) 
	if option=='full_boss':
		class_best = data['DR16_CLASS'] 
	else:
		class_best = data['merged_class']
	print(np.unique(class_best, return_counts=True))
	#
	high_conf = (data['RXS_ExiML'] >= EXI_ML_min )
	targeted  = (high_conf) & (data['SDSS_FIBER2MAG_i']>=MIN_SDSS_FIBER2MAG_i) & (data['SDSS_FIBER2MAG_i']<=MAX_SDSS_FIBER2MAG_i) & (data['SDSS_MODELMAG_i']>=MIN_SDSS_MODELMAG_i)    
	observed  = (targeted) & (data['DR16_MEMBER']==1)
	c1 = (observed) & ((data['Z_BEST']>-0.5) | ((data['DR16_Z']>-0.5) & (data['DR16_Z_ERR']>0)))
	c2 = (c1) & (data['CONF_BEST']==3)
	c3 = (c1) & (data['CONF_BEST']==2) & ((class_best=='BLAZAR')|(class_best=='BLLAC'))
	c4 = (c1) & (data['DR16_SN_MEDIAN_ALL']>=2) & (data['DR16_ZWARNING']==0 )
	c5 = (c1) & (data['CONF_BEST']==2) & (data['DR16_ZWARNING']==0 ) & (data['VI_REINSPECT_FLAG'] == 0) & (data['VI_NINSPECTORS']>2)
	c6 = (c1) & (data['CONF_BEST']==2) & (data['DR16_ZWARNING']==0 ) & (data['VI_AM_RECONCILED']==1)
	idZ =  (c2) | (c3) | (c4) | (c5) | (c6) 
	blazars_noZ = (idZ) & ( (class_best=='BLAZAR') | (class_best=='BLLAC')) & (data['CONF_BEST']<3)
	goodZ = (idZ) & (blazars_noZ == False)
	redmapper_cluster = (goodZ) & (data['REDMAPPER_Separation']<60) & (abs(data['DR16_Z'] - data['REDMAPPER_Z_LAMBDA'])<0.01)
	spiders_cluster = (goodZ) & (abs(data['SPIDERSCODEX_SCREEN_CLUZSPEC']-data['Z_BEST'])<0.01)
	clusters = (redmapper_cluster) | (spiders_cluster)
	stars = (class_best=='STAR')
	agnZ = (goodZ) & (clusters==False) &  (stars==False)
	bl  = (agnZ) & ( (class_best=='BLAGN')  | (class_best=='QSO') )
	other = (agnZ) & ( (class_best=="BLLAC") | ( class_best== "BALQSO")  | ( class_best== "BLAZAR")  | ( class_best== "NONE")  | ( class_best== "QSO_BAL") )  
	nl  = (agnZ) & ( bl==False) & (other==False)
	print(len(data), len(targeted), len(observed), len(goodZ), len(idZ), len(agnZ), len(clusters), len(blazars_noZ), len(stars), len(bl), len(nl) )
	return data, high_conf, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, bl, nl, class_best

class Catalog: pass

data, high_conf, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, bl, nl, class_best = get_arrays_xmmsl2(path_2_cat_xmmsl = path_2_cat_xmmsl, EXI_ML_min = EXI_ML_min)

data_XMM = Catalog()
data_XMM.data = data
data_XMM.high_conf = high_conf
data_XMM.targeted = targeted
data_XMM.observed = observed
data_XMM.goodZ = goodZ
data_XMM.idZ = idZ
data_XMM.agnZ = agnZ
data_XMM.clusters = clusters
data_XMM.blazars_noZ = blazars_noZ
data_XMM.stars = stars
data_XMM.bl = bl
data_XMM.nl = nl
data_XMM.class_best = class_best

if write_table:
	new_tab = Table(data)
	new_tab.add_column(Column(name='is_target', data=targeted ))
	new_tab.add_column(Column(name='is_observed', data=observed ))
	new_tab.add_column(Column(name='is_goodZ', data=goodZ    ))
	new_tab.add_column(Column(name='is_identified', data=idZ      ))
	new_tab.add_column(Column(name='is_agn', data=agnZ     ))
	new_tab.add_column(Column(name='is_cluster', data=clusters ))
	new_tab.add_column(Column(name='is_blazar', data=blazars_noZ ))
	new_tab.add_column(Column(name='is_star', data=stars ))
	new_tab.add_column(Column(name='is_bl', data=bl ))
	new_tab.add_column(Column(name='is_nl', data=nl ))
	new_tab.add_column(Column(name='merged_class', data=class_best ))
	path_2_cat_xmmsl_out = os.path.join(catalog_dir, 'XMMSL2',  'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_REDMAPPER_SPIDERSCODEX_XID_classifications.fits')
	new_tab.write(path_2_cat_xmmsl_out)

if update_table:
	print('updates')
	new_tab = Table(data)
	new_tab['is_bl'] = bl
	new_tab['is_nl'] = nl
	print( np.unique( new_tab['merged_class'][new_tab['is_bl']], return_counts = True))
	print( np.unique( new_tab['merged_class'][new_tab['is_nl']], return_counts = True))
	new_tab.write(path_2_cat_xmmsl, overwrite=True)
	

data, high_conf, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, bl, nl, class_best = get_arrays_2rxs(path_2_cat_2RXS = path_2_cat_2RXS, EXI_ML_min = EXI_ML_min)

data_RXS = Catalog()
data_RXS.data = data
data_RXS.high_conf = high_conf
data_RXS.targeted = targeted
data_RXS.observed = observed
data_RXS.goodZ    = goodZ
data_RXS.idZ      = idZ
data_RXS.agnZ     = agnZ
data_RXS.clusters = clusters
data_RXS.blazars_noZ = blazars_noZ
data_RXS.stars = stars
data_RXS.bl = bl
data_RXS.nl = nl
data_RXS.class_best = class_best

if write_table:
	new_tab = Table(data)
	new_tab.add_column(Column(name='is_target', data=targeted ))
	new_tab.add_column(Column(name='is_observed', data=observed ))
	new_tab.add_column(Column(name='is_goodZ', data=goodZ    ))
	new_tab.add_column(Column(name='is_identified', data=idZ      ))
	new_tab.add_column(Column(name='is_agn', data=agnZ     ))
	new_tab.add_column(Column(name='is_cluster', data=clusters ))
	new_tab.add_column(Column(name='is_blazar', data=blazars_noZ ))
	new_tab.add_column(Column(name='is_star', data=stars ))
	new_tab.add_column(Column(name='is_bl', data=bl ))
	new_tab.add_column(Column(name='is_nl', data=nl ))
	new_tab.add_column(Column(name='merged_class', data=class_best ))
	path_2_cat_2RXS_out = os.path.join(catalog_dir, '2RXS', 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER_SPIDERSCODEX_XID_classifications.fits')
	new_tab.write(path_2_cat_2RXS_out)

if update_table:
	print('updates')
	new_tab = Table(data)
	new_tab['is_bl'] = bl
	new_tab['is_nl'] = nl
	print( np.unique( new_tab['merged_class'][new_tab['is_bl']], return_counts = True))
	print( np.unique( new_tab['merged_class'][new_tab['is_nl']], return_counts = True))
	new_tab.write(path_2_cat_2RXS, overwrite=True)
	
