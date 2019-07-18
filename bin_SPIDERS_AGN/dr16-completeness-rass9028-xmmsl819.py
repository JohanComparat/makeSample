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

catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )
figure_dir = os.path.join(catalog_dir, 'figures' )

if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)

'RASS_3XMM_COLOURMAG_PIX_smoothed2pix_bestWISE_withSDSSv5bbest.fits'
'XMMSLEW_3XMM_COLOURMAG_PIX_smoothed2pix_bestWISE_withSDSSv5bbest.fits'

path_2_rass_cat = os.path.join(catalog_dir, 'spiderstargetAGN-SPIDERS_RASS_AGN-v2.1_v5_13_0_sdss_26_VERON_2QZ_MASKED.fits')
path_2_xmmsl_cat = os.path.join(catalog_dir, 'spiderstargetAGN-SPIDERS_XMMSL_AGN-v3.1_v5_13_0_sdss_26_VERON_2QZ_MASKED.fits')

hd_rass = fits.open(path_2_rass_cat)[1].data
hd_xmmsl = fits.open(path_2_xmmsl_cat)[1].data

data = hd_rass

print( len(data) )
print( len(data[(data['in_BOSS_v5_13_0'])]), len(data[(data['in_BOSS_v5_13_0'])])*1./len(data) )
print( len(data[(data['in_BOSS_v5_13_0'])&(data['ZWARNING']==0)]), len(data[(data['in_BOSS_v5_13_0'])&(data['ZWARNING']==0)])*1./len(data)  )
observed = (data['in_BOSS_v5_13_0'])
goodZ = (data['in_BOSS_v5_13_0'])&(data['ZWARNING']==0)
# fraction of target observed and fractino of successful redshifts as a function of quantities. 

qty = 'XRAY_FLUX'
YY = np.log10(data[qty])

bins = np.arange(np.min(YY), np.max(YY), (np.max(YY)-np.min(YY))/20.)
x_bins = (bins[1:]+bins[:-1])/2.

out_all = np.histogram( YY, bins=bins )[0]
out_obs = np.histogram( YY[observed], bins=bins )[0]
out_zzz = np.histogram( YY[goodZ], bins=bins )[0]


fig_out = os.path.join(figure_dir, qty+'_hist.png')

p.figure(1, (5.5,5.5))
p.axes([0.17, 0.17, 0.78, 0.75])
p.tight_layout()
p.plot(x_bins, out_all/np.sum(out_all), label='pdf')
#
mid_line = out_zzz*1./out_obs
err_line = out_zzz**(-0.5)
p.fill_between(x_bins, y1=mid_line*(1-err_line), y2=mid_line*(1+err_line), label='good Z/observed', alpha=0.5)
#
mid_line = out_obs*1./out_all
err_line = out_obs**(-0.5)
p.fill_between(x_bins, y1=mid_line*(1-err_line), y2=mid_line*(1+err_line), label='observed/all', alpha=0.5)
#
p.xlabel(qty)#'Observed fraction per '+pixel_area_str+r' deg$^2$ pixels')
p.ylabel('fraction')
p.grid()
p.ylim((0,1))
#p.yscale('log')
p.legend(frameon=False, loc=4)
p.savefig(fig_out)
p.clf()


# NSIDE=2^5=32
# 
nside_int = int(sys.argv[1])
nside = 2**nside_int
nside_str = 'hp'+str(nside_int)
pixel_area = healpy.nside2pixarea(nside, degrees=True)
pixel_area_str = str(np.round(pixel_area,2))
uniq_all = np.unique(data[nside_str], return_counts=True)
uniq_observed = np.unique(data[nside_str][data['in_BOSS_v5_13_0']], return_counts=True)
uniq_goodZ = np.unique(data[nside_str][(data['in_BOSS_v5_13_0'])&(data['ZWARNING']==0)], return_counts=True)

N_pixels_with_targets = len(uniq_all[0])
N_pixels_with_observed_targets = len(uniq_observed[0])
print('N_pixels_with_targets', N_pixels_with_targets)
print('N_pixels_with_observed_targets', N_pixels_with_observed_targets)
N_pixels_with_goodZ = len(uniq_goodZ[0])
print('N_pixels_with_goodZ', N_pixels_with_goodZ)

all_targets = interp1d(uniq_all[0], uniq_all[1]*1.)
completeness = np.round( uniq_observed[1]/all_targets(uniq_observed[0]), 3)
success_rate = np.round( uniq_goodZ[1]/all_targets(uniq_goodZ[0]), 3)

N_pixels_with_100pc_observed_targets = len(completeness[completeness==1])
print( 'N_pixels_with_100pc_observed_targets', N_pixels_with_100pc_observed_targets )
print( 'area per pixel ', pixel_area, 'deg2')
print( 'area ', N_pixels_with_100pc_observed_targets*pixel_area, 'deg2')

fig_out = os.path.join(figure_dir, 'completeness_hist_'+nside_str+'.png')

p.figure(1, (5.5,5.5))
p.axes([0.17, 0.17, 0.78, 0.75])
p.tight_layout()
#p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2)
out = p.hist(completeness, bins = np.arange(0.,0.99, 0.01, ), histtype='step', rasterized=True, lw=3)
#p.title(baseName)
p.xlabel('Observed fraction per '+pixel_area_str+r' deg$^2$ pixels')
p.ylabel('Number of pixels')
p.grid()
#p.yscale('log')
#p.legend(frameon=False, loc=0)
p.savefig(fig_out)
p.clf()

fig_out = os.path.join(figure_dir, 'successZ_hist_'+nside_str+'.png')

p.figure(1, (5.5,5.5))
p.axes([0.17, 0.17, 0.78, 0.75])
p.tight_layout()
#p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2)
out = p.hist(success_rate, bins = np.arange(0.,0.99, 0.01, ), histtype='step', rasterized=True, lw=3)
#p.title(baseName)
p.xlabel('Observed fraction per '+pixel_area_str+r' deg$^2$ pixels')
p.ylabel('Number of pixels')
p.grid()
#p.yscale('log')
#p.legend(frameon=False, loc=0)
p.savefig(fig_out)
p.clf()

NN = out[0].astype('int')
print('sum(NN)', np.sum(NN))
print('histogram', NN)
print('histogram fraction', np.round(100*NN*1./N_pixels_with_observed_targets,1))

not_observed = np.in1d(uniq_all[0], uniq_observed[0], invert=True)

ra_all, dec_all = healpy.pix2ang(nside, uniq_all[0][not_observed], nest=True, lonlat=True)
ra, dec = healpy.pix2ang(nside, uniq_observed[0], nest=True, lonlat=True)

fig_out = os.path.join(figure_dir, 'completeness_ra_dec_'+nside_str+'.png')

p.figure(1, (6.5,5.5))
p.axes([0.17, 0.17, 0.78, 0.75])
p.tight_layout()
#p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2)
p.plot(ra_all, dec_all, 'k,', label='not observed')
p.scatter(ra, dec, c=completeness, s=16-2*nside_int, rasterized=True)
p.colorbar(shrink=0.9)
p.xlabel('R.A. [degree]')
p.ylabel('Dec. [degree]')
p.grid()
p.title('Observed fraction per '+pixel_area_str+r' deg$^2$ pixels')
#p.yscale('log')
#p.legend(frameon=False, loc=0)
p.savefig(fig_out)
p.clf()

ra, dec = healpy.pix2ang(nside, uniq_goodZ[0], nest=True, lonlat=True)

fig_out = os.path.join(figure_dir, 'successZ_ra_dec_'+nside_str+'.png')

p.figure(1, (6.5,5.5))
p.axes([0.17, 0.17, 0.78, 0.75])
p.tight_layout()
#p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2)
p.scatter(ra, dec, c=success_rate, s=16-2*nside_int, rasterized=True)
p.colorbar(shrink=0.9)
p.xlabel('R.A. [degree]')
p.ylabel('Dec. [degree]')
p.grid()
p.title('Success rate per '+pixel_area_str+r' deg$^2$ pixels')
#p.yscale('log')
#p.legend(frameon=False, loc=0)
p.savefig(fig_out)
p.clf()

