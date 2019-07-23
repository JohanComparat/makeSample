"""

python dr16-completeness-VAC-2RXS.py 4 6p5
python dr16-completeness-VAC-2RXS.py 4 10

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

nside_int = int(sys.argv[1])
nside = 2**nside_int
nside_str = 'hp'+str(nside_int)
pixel_area = healpy.nside2pixarea(nside, degrees=True)
pixel_area_str = str(np.round(pixel_area,2))

catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )
figure_dir = os.path.join(catalog_dir, 'figures_VAC' )

if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)

#path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_2RXS_DR16.fits')

path_2_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits')
# contains 135,259
#hd0=fits.open(os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits'))
# contains 21,288
path_2_cat_dr16 = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')

MIN_SDSS_FIBER2MAG_i = 17.0
MAX_SDSS_FIBER2MAG_i = 22.5
MIN_SDSS_MODELMAG_i  = 16.0

EXI_ML_min_str = sys.argv[2]

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
	return data[(high_conf)], ra[(high_conf)], dec[(high_conf)], targeted[(high_conf)], observed[(high_conf)], goodZ[(high_conf)]

#data['RXS_ExtML'    ]
#data['RXS_ExiML1RXS'] 
#data['RXS_ExiML2RXS']

# 2RXS
qty = 'RXS_SRC_FLUX'
# XMMSL2
#qty = 'XMMSL2_FLUX_B6'

fig_out = os.path.join(figure_dir, prefix+qty+'_hist_dr16.png')


XMIN = -13. 
XMAX = -10
DX=0.5
bins = np.arange(XMIN, XMAX+DX, DX)# (np.max(YY)-np.min(YY))/10.)
x_bins = (bins[1:]+bins[:-1])/2.

p.figure(1, (5.5,5.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()

# DR16 footprint
print('======================================')
print('======================================')
print('DR16 footprint')
print('======================================')
print('======================================')
data, ra, dec, targeted, observed, goodZ = get_arrays(path_2_cat_dr16)
print( len(data) )
print( len(data[targeted]), len(data[targeted])*1./len(data) )
print( len(data[observed]), len(data[observed])*1./len(data[targeted]) )
print( len(data[goodZ]), len(data[goodZ])*1./len(data[observed]) )
# fraction of target observed and fractino of successful redshifts as a function of quantities. 
YY = np.log10(data[qty])

out_all = np.histogram( YY, bins=bins )[0]
out_tar = np.histogram( YY[targeted], bins=bins )[0]
out_obs = np.histogram( YY[observed], bins=bins )[0]
out_zzz = np.histogram( YY[goodZ], bins=bins )[0]

p.plot(x_bins, out_all/np.sum(out_all), label='pdf DR16')
#
mid_line = out_tar*1./out_all
err_line = out_tar**(-0.5)
p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='T/A', fmt='none', lw=3)
#
mid_line = out_obs*1./out_tar
err_line = out_obs**(-0.5)
p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='O/T', fmt='none', lw=2)
#
mid_line = out_zzz*1./out_obs
err_line = out_zzz**(-0.5)
p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='Z/O', fmt='none')

p.xlabel(r'$\log_{10}(F_X [erg.\, cm^{-2}. s^{-1}])$')
p.ylabel('fraction')
p.grid()
p.xlim((XMIN, XMAX))
p.ylim((0,1))
#p.yscale('log')
#p.legend(frameon=True, loc=7)
p.savefig(fig_out)
p.clf()


fig_out = os.path.join(figure_dir, prefix+qty+'_hist.png')

p.figure(1, (5.5,5.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()

# FULL footprint
print('======================================')
print('======================================')
print('FULL footprint')
print('======================================')
print('======================================')
data, ra, dec, targeted, observed, goodZ = get_arrays(path_2_cat)
print( 'A', len(data) )
print( 'T',len(data[targeted]), len(data[targeted])*1./len(data) )
print( 'O',len(data[observed]), len(data[observed])*1./len(data[targeted]) )
print( 'Z',len(data[goodZ]), len(data[goodZ])*1./len(data[observed]) )
# fraction of target observed and fractino of successful redshifts as a function of quantities. 
YY = np.log10(data[qty])

out_all = np.histogram( YY, bins=bins )[0]
out_tar = np.histogram( YY[targeted], bins=bins )[0]
out_obs = np.histogram( YY[observed], bins=bins )[0]
out_zzz = np.histogram( YY[goodZ], bins=bins )[0]


p.plot(x_bins, out_all/np.sum(out_all), label='pdf')
#
mid_line = out_tar*1./out_all
err_line = out_tar**(-0.5)
p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='T/A', fmt='none', lw=3)
#
mid_line = out_obs*1./out_tar
err_line = out_obs**(-0.5)
p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='O/T', fmt='none', lw=2)
#
mid_line = out_zzz*1./out_obs
err_line = out_zzz**(-0.5)
p.errorbar(x_bins, y = mid_line, xerr=DX/2., yerr=err_line, label='Z/O', fmt='none')

p.xlabel(qty)#'Observed fraction per '+pixel_area_str+r' deg$^2$ pixels')
p.xlabel(r'$\log_{10}(F_X [erg.\, cm^{-2}. s^{-1}])$')
p.ylabel('fraction')
p.grid()
p.xlim((XMIN, XMAX))
p.ylim((0,1))
#p.yscale('log')
p.legend(frameon=True, loc=3)
p.savefig(fig_out)
p.clf()


# NSIDE=2^5=32
# 

hpx = healpy.ang2pix(nside, ra, dec, nest=True, lonlat=True)
uniq_all = np.unique(hpx, return_counts=True)
uniq_targeted = np.unique(hpx[targeted], return_counts=True)
uniq_observed = np.unique(hpx[observed], return_counts=True)
uniq_goodZ = np.unique(hpx[goodZ], return_counts=True)

N_pixels_with_targets = len(uniq_all[0])
print('N_pixels_with_targets', N_pixels_with_targets)
N_pixels_with_observed_targets = len(uniq_observed[0])
print('N_pixels_with_observed_targets', N_pixels_with_observed_targets)
N_pixels_with_goodZ = len(uniq_goodZ[0])
print('N_pixels_with_goodZ', N_pixels_with_goodZ)

all_targets = interp1d(uniq_targeted[0], uniq_targeted[1]*1.)
completeness = np.round( uniq_observed[1]/all_targets(uniq_observed[0]), 3)
success_rate = np.round( uniq_goodZ[1]/all_targets(uniq_goodZ[0]), 3)

fig_out = os.path.join(figure_dir, prefix+'successZ_hist_'+nside_str+'.png')

p.figure(1, (5.5,5.5))
p.axes([0.17, 0.17, 0.78, 0.75])
p.tight_layout()
out1 = p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2, label='O/T')
out  = p.hist(success_rate, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=3, label='Z/O')
#p.title(baseName)
p.xlabel('Observed fraction per '+pixel_area_str+r' deg$^2$ pixels')
p.ylabel('Number of pixels')
p.grid()
#p.yscale('log')
p.legend(frameon=False, loc=0)
p.savefig(fig_out)
p.clf()

not_observed = np.in1d(uniq_all[0], uniq_observed[0], invert=True)

ra_all, dec_all = healpy.pix2ang(nside, uniq_all[0][not_observed], nest=True, lonlat=True)
ra_hp, dec_hp = healpy.pix2ang(nside, uniq_observed[0], nest=True, lonlat=True)

p2map = os.path.join(catalog_dir, 'DR16_footprint_hpx10.fits')
p2map = os.path.join(catalog_dir, 'DR16_wSEQUELS_footprint_clipped_complete_hpx10.fits')
p2map = os.path.join(catalog_dir, 'DR16_wSEQUELS_footprint_hpx10.fits')
hd=fits.open(p2map)

NSIDE = 1024
all_pix = np.arange(healpy.nside2npix(NSIDE))
ra_mask, dec_mask = healpy.pix2ang(NSIDE, all_pix[hd[1].data['MASK'].astype('bool')], nest=True, lonlat=True)

fig_out = os.path.join(figure_dir, prefix+'successZ_ra_dec_'+nside_str+'.png')

p.figure(1, (7.,5.5))
p.axes([0.17, 0.17, 0.78, 0.75])
p.tight_layout()
#p.plot(ra_mask, dec_mask, 'ko', rasterized=True, alpha=0.0015)
#p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2)
p.scatter(ra_hp, dec_hp, c=completeness, s=20, marker='o', rasterized=True, cmap='RdBu', label='exiML>'+str(np.round(EXI_ML_min,1)))#, alpha=0.75)
p.colorbar(shrink=0.9)
p.xlabel('R.A. [degree]')
p.ylabel('Dec. [degree]')
p.grid()
p.title('Fraction observed per '+pixel_area_str+r' deg$^2$ pixels')
#p.yscale('log')
p.legend(frameon=False, loc=0)
p.savefig(fig_out)
p.clf()

fig_out = os.path.join(figure_dir, 'eBOSS-dr16-wSequels.png')

p.figure(1, (6.5,5.5))
p.axes([0.17, 0.17, 0.78, 0.75])
p.tight_layout()
p.plot(ra_mask, dec_mask, 'ko', rasterized=True)
p.xlabel('R.A. [degree]')
p.ylabel('Dec. [degree]')
p.grid()
#p.yscale('log')
#p.legend(frameon=False, loc=0)
p.savefig(fig_out)
p.clf()


