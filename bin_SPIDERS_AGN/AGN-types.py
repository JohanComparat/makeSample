"""

python AGN-types.py

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

from sklearn.neighbors import BallTree, DistanceMetric
from astropy.table import Table,unique,Column
from math import radians, cos, sin, asin, sqrt, pi

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

agn_clustering_dir = '/data36s/comparat/AGN_clustering'
agn_clustering_dir = '/home/comparat/data/AGN_clustering'

#nside_int = int(sys.argv[1])
#nside = 2**nside_int
#nside_str = 'hp'+str(nside_int)
#pixel_area = healpy.nside2pixarea(nside, degrees=True)
#pixel_area_str = str(np.round(pixel_area,2))

catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )
figure_dir = os.path.join(catalog_dir, 'figures_VAC' )

if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)


path_2_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits')
path_2_cat_dr16 = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')

MIN_SDSS_FIBER2MAG_i = 17.0
MAX_SDSS_FIBER2MAG_i = 22.5
MIN_SDSS_MODELMAG_i  = 16.0

EXI_ML_min_str = '6p5' # sys.argv[2]

if EXI_ML_min_str=='6p5':
	EXI_ML_min=6.5
if EXI_ML_min_str=='10':
	EXI_ML_min=10

prefix='2RXS_ExiML'+EXI_ML_min_str+'_'

#prefix='XMMSL2_ExiML10_'
#prefix='XMMSL2_ExiML6p5_'

#def get_arrays(path_2_cat, EXI_ML_min = 10.):
def get_arrays(path_2_cat, EXI_ML_min = EXI_ML_min):
	hd_rass_i = fits.open(path_2_cat)[1].data
	targets = (hd_rass_i['NWAY_match_flag']!=2) & (hd_rass_i['RXS_IN_BOSS']==1) & (hd_rass_i['FLAG_SDSSv5b_best']==1) & (hd_rass_i['NWAY_p_any']>=0.01) & (hd_rass_i['NUM_SDSS']>=1)
	#  NWAY_match_flag!=2 &&RXS_IN_BOSS==1 &&FLAG_SDSSv5b_best==1 &&
	hd_rass = fits.open(path_2_cat)[1].data[targets]
	#print(np.unique(hd_rass['CLASS_BEST']))
	data = hd_rass
	ra = data['RXS_RAJ2000']
	dec = data['RXS_DEJ2000']
	high_conf = (data['RXS_ExiML'] >= EXI_ML_min )
	targeted  = (high_conf) & (data['SDSS_FIBER2MAG_i']>=MIN_SDSS_FIBER2MAG_i) & (data['SDSS_FIBER2MAG_i']<=MAX_SDSS_FIBER2MAG_i) & (data['SDSS_MODELMAG_i']>=MIN_SDSS_MODELMAG_i)    
	observed  = (targeted) & (data['DR16_MEMBER']==1)
	goodZ     = (observed) & (data['CONF_BEST']==3)
	print(len(ra), len(ra[(high_conf)]))
	return data[observed], goodZ[observed]

data, goodZ = get_arrays(path_2_cat_dr16)

all_z = (data['CONF_BEST']==3)|(data['CONF_BEST']==2)

# We distinguish four types of AGN: broad line, narrow line, broad absorption line, blazar
#if CLASS\_BEST = BLAGN and DR16\_CLASS = GALAXY, then group with GALAXY
#if CLASS\_BEST = BLAGN and DR16\_CLASS = QSO, then group with QSO


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

qty='SDSS_MODELMAG_i'

fig_out = os.path.join(figure_dir, prefix+'_'+qty+'_types_dr16.png')

#YY = data[qty][good]#'DR16Q_SN_MEDIAN_ALL']

XMIN = 16 
XMAX = 30
DX=1.0 #0.5
bins = np.arange(XMIN, XMAX+DX, DX)# (np.max(YY)-np.min(YY))/10.)
x_bins = (bins[1:]+bins[:-1])/2.

p.figure(1, (5.5,5.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()

out_bl  = np.histogram( data[qty][bl ], bins=bins)[0]
out_nl  = np.histogram( data[qty][nl ], bins=bins)[0]
out_ar  = np.histogram( data[qty][ar ], bins=bins)[0]
out_bal = np.histogram( data[qty][bal], bins=bins)[0]
out_no  = np.histogram( data[qty][no ], bins=bins)[0]
out_st  = np.histogram( data[qty][st ], bins=bins)[0]
#out_ = np.histogram( data[qty][], bins=bins)[0]

#p.plot(x_bins, out_all/np.sum(out_all),lw=3, label='pdf')
#p.plot(x_bins, out_cb0/np.sum(out_cb0),lw=2, label='0')
p.plot(x_bins, out_bl /np.sum(out_bl ),lw=3, label='bl , N='+str(int(np.sum(out_bl ))) )
p.plot(x_bins, out_nl /np.sum(out_nl ),lw=3, label='nl , N='+str(int(np.sum(out_nl ))) )
p.plot(x_bins, out_ar /np.sum(out_ar ),lw=3, label='blazar , N='+str(int(np.sum(out_ar ))), ls='dashed' )
p.plot(x_bins, out_bal/np.sum(out_bal),lw=3, label='bal, N='+str(int(np.sum(out_bal))), ls='dashed' )
p.plot(x_bins, out_no /np.sum(out_no ),lw=3, label='none , N='+str(int(np.sum(out_no ))), ls='dotted' )
p.plot(x_bins, out_st /np.sum(out_st ),lw=3, label='star , N='+str(int(np.sum(out_st ))), ls='dotted' )
#
#mid_line = out_tar*1./out_all
#err_line = out_tar**(-0.5)
#p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='T/A', fmt='none', lw=3)
##
#mid_line = out_obs*1./out_tar
#err_line = out_obs**(-0.5)
#p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='O/T', fmt='none', lw=2)
##
#mid_line = out_zzz*1./out_obs
#err_line = out_zzz**(-0.5)
#p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='Z/O', fmt='none')

p.xlabel(qty)
p.ylabel('distribution')
p.grid()
p.xlim((bins.min(), 23 )) # np.max(x_bins[out_all>1]) + DX ))
#p.ylim((0,1))
#p.xscale('log')
p.legend(frameon=True, loc=0)
p.savefig(fig_out)
p.clf()


qty='SDSS_FIBER2MAG_i'

fig_out = os.path.join(figure_dir, prefix+'_'+qty+'_types_dr16.png')

#YY = data[qty][good]#'DR16Q_SN_MEDIAN_ALL']

XMIN = 16 
XMAX = 30
DX=1.0 #0.5
bins = np.arange(XMIN, XMAX+DX, DX)# (np.max(YY)-np.min(YY))/10.)
x_bins = (bins[1:]+bins[:-1])/2.

p.figure(1, (5.5,5.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()

out_bl  = np.histogram( data[qty][bl ], bins=bins)[0]
out_nl  = np.histogram( data[qty][nl ], bins=bins)[0]
out_ar  = np.histogram( data[qty][ar ], bins=bins)[0]
out_bal = np.histogram( data[qty][bal], bins=bins)[0]
out_no  = np.histogram( data[qty][no ], bins=bins)[0]
out_st  = np.histogram( data[qty][st ], bins=bins)[0]
#out_ = np.histogram( data[qty][], bins=bins)[0]

#p.plot(x_bins, out_all/np.sum(out_all),lw=3, label='pdf')
#p.plot(x_bins, out_cb0/np.sum(out_cb0),lw=2, label='0')
p.plot(x_bins, out_bl /np.sum(out_bl ),lw=3, label='bl , N='+str(int(np.sum(out_bl ))) )
p.plot(x_bins, out_nl /np.sum(out_nl ),lw=3, label='nl , N='+str(int(np.sum(out_nl ))) )
p.plot(x_bins, out_ar /np.sum(out_ar ),lw=3, label='blazar , N='+str(int(np.sum(out_ar ))), ls='dashed' )
p.plot(x_bins, out_bal/np.sum(out_bal),lw=3, label='bal, N='+str(int(np.sum(out_bal))), ls='dashed' )
p.plot(x_bins, out_no /np.sum(out_no ),lw=3, label='none , N='+str(int(np.sum(out_no ))), ls='dotted' )
p.plot(x_bins, out_st /np.sum(out_st ),lw=3, label='star , N='+str(int(np.sum(out_st ))), ls='dotted' )
#
#mid_line = out_tar*1./out_all
#err_line = out_tar**(-0.5)
#p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='T/A', fmt='none', lw=3)
##
#mid_line = out_obs*1./out_tar
#err_line = out_obs**(-0.5)
#p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='O/T', fmt='none', lw=2)
##
#mid_line = out_zzz*1./out_obs
#err_line = out_zzz**(-0.5)
#p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='Z/O', fmt='none')

p.xlabel(qty)
p.ylabel('distribution')
p.grid()
p.xlim((bins.min(), 23 )) #np.max(x_bins[out_all>1]) + DX ))
#p.ylim((0,1))
#p.xscale('log')
p.legend(frameon=True, loc=0)
p.savefig(fig_out)
p.clf()


qty='Z_BEST'

fig_out = os.path.join(figure_dir, prefix+'_'+qty+'_types_dr16.png')

#YY = data[qty][good]#'DR16Q_SN_MEDIAN_ALL']

XMIN = 0. 
XMAX = 4.
DX= 0.25 #0.5
bins = np.arange(XMIN, XMAX+DX, DX)# (np.max(YY)-np.min(YY))/10.)
x_bins = (bins[1:]+bins[:-1])/2.

p.figure(1, (5.5,5.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()

out_bl  = np.histogram( data[qty][bl ], bins=bins)[0]
out_nl  = np.histogram( data[qty][nl ], bins=bins)[0]
out_ar  = np.histogram( data[qty][ar ], bins=bins)[0]
out_bal = np.histogram( data[qty][bal], bins=bins)[0]
out_no  = np.histogram( data[qty][no ], bins=bins)[0]
out_st  = np.histogram( data[qty][st ], bins=bins)[0]
#out_ = np.histogram( data[qty][], bins=bins)[0]

#p.plot(x_bins, out_all/np.sum(out_all),lw=3, label='pdf')
#p.plot(x_bins, out_cb0/np.sum(out_cb0),lw=2, label='0')
p.hist(data[qty][bl ], bins=bins, histtype='step',lw=3, label='bl , N='+str(int(np.sum(out_bl ))) )
p.hist(data[qty][nl ], bins=bins, histtype='step',lw=3, label='nl , N='+str(int(np.sum(out_nl ))) )
p.hist(data[qty][ar ], bins=bins, histtype='step',lw=3, label='blazar , N='+str(int(np.sum(out_ar ))), ls='dashed' )
p.hist(data[qty][bal], bins=bins, histtype='step',lw=3, label='bal, N='+str(int(np.sum(out_bal))), ls='dashed' )
p.hist(data[qty][no ], bins=bins, histtype='step',lw=4, label='none , N='+str(int(np.sum(out_no ))), ls='dotted' )
p.hist(data[qty][st ], bins=bins, histtype='step',lw=4, label='star , N='+str(int(np.sum(out_st ))), ls='dotted' )
#
#mid_line = out_tar*1./out_all
#err_line = out_tar**(-0.5)
#p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='T/A', fmt='none', lw=3)
##
#mid_line = out_obs*1./out_tar
#err_line = out_obs**(-0.5)
#p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='O/T', fmt='none', lw=2)
##
#mid_line = out_zzz*1./out_obs
#err_line = out_zzz**(-0.5)
#p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='Z/O', fmt='none')

p.xlabel('redshift')
p.ylabel('distribution')
p.grid()
#p.xlim((bins.min(), 23 )) #np.max(x_bins[out_all>1]) + DX ))
#p.ylim((0,1))
p.yscale('log')
p.legend(frameon=True, loc=0)
p.savefig(fig_out)
p.clf()
