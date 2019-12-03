
from catalog_lib import *

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

nside_int = int(sys.argv[1])
nside = 2**nside_int
nside_str = 'hp'+str(nside_int)
pixel_area = healpy.nside2pixarea(nside, degrees=True)
pixel_area_str = str(np.round(pixel_area,2))

EXI_ML_min_str = '10'
EXI_ML_min = 10

prefix='XMMSL2_ExiML'+EXI_ML_min_str+'_'

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


qty = 'XMMSL2_FLUX_B6'

fig_out = os.path.join(figure_dir, prefix+qty+'_hist_dr16.png')


XMIN = -13. 
XMAX = -10
DX=1.
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


ExiML = np.transpose([
	data_XMM.data['XMMSL2_DET_ML_B6'] , 
	data_XMM.data['XMMSL2_DET_ML_B7'] , 
	data_XMM.data['XMMSL2_DET_ML_B8'] ])
ExiML[np.isnan(ExiML)]=0
ExiML_arr = np.max(ExiML,axis=1) 

print('All         &', len(data_XMM.data) )
print('Targets     &', len(data_XMM.data[data_XMM.targeted]) )#, len(data_XMM.data[data_XMM.targeted])*1./len(data_XMM.data) )
print('Observed    &', len(data_XMM.data[data_XMM.observed]) )#, len(data_XMM.data[data_XMM.observed])*1./len(data_XMM.data[data_XMM.targeted]) )
print('Identified  &', len(data_XMM.data[data_XMM.idZ])      )#, len(data_XMM.data[data_XMM.idZ])*1./len(data_XMM.data[data_XMM.observed]) )
print('Blazar no Z &', len(data_XMM.data[data_XMM.blazars_noZ]))#, len(data_XMM.data[data_XMM.blazars_noZ])*1./len(data_XMM.data[data_XMM.observed]) )
print('good Z      &', len(data_XMM.data[data_XMM.goodZ])    )#, len(data_XMM.data[data_XMM.goodZ])*1./len(data_XMM.data[data_XMM.observed]) )
print('AGN         &', len(data_XMM.data[data_XMM.agnZ])     )#, len(data_XMM.data[data_XMM.agnZ])*1./len(data_XMM.data[data_XMM.observed]) )
print('Cluster     &', len(data_XMM.data[data_XMM.clusters]) )#, len(data_XMM.data[data_XMM.clusters])*1./len(data_XMM.data[data_XMM.observed]) )
print('Star        &', len(data_XMM.data[data_XMM.stars])    )#, len(data_XMM.data[data_XMM.stars])*1./len(data_XMM.data[data_XMM.observed]) )


dz = 0.1
z_min=0.0005 
redshift = data_XMM.data['DR16_Z']
z_max = np.max(redshift[data_XMM.agnZ]) # np.max([np.max(redshift),np.max(Z_XMMSL2)])
zs = np.arange(z_min, z_max + dz, dz)
xzs = (zs[1:]+zs[:-1])/2.

NN_xmmsl2 = np.histogram(redshift[data_XMM.agnZ], bins=zs)[0]
#density = NN_xmmsl2/area_eboss_dr16#/np.max(NN)

NN_clusters = np.histogram(redshift[data_XMM.clusters], bins=zs)[0]
#density = NN_clusters/area_eboss_dr16#/np.max(NN)

np.savetxt(os.path.join( figure_dir, "histogram_redshift.txt"), np.transpose([zs[:-1], zs[1:], NN_xmmsl2]), delimiter=" & ", newline=" \\\\ \n", fmt='%10.1f')

out = np.unique(data_XMM.data['CLASS_BEST'][data_XMM.clusters], return_counts=True)
print(out)

agn_in_cluster = (data_XMM.clusters) & (data_XMM.data['CLASS_BEST']!="GALAXY")
for el in np.transpose([data_XMM.data['ALLW_RA'][agn_in_cluster],data_XMM.data['ALLW_DEC'][agn_in_cluster], data_XMM.data['DR16_Z'][agn_in_cluster], data_XMM.data['REDMAPPER_Z_LAMBDA'][agn_in_cluster], data_XMM.data['REDMAPPER_Separation'][agn_in_cluster]]):
	print(np.round(el,5))


#sys.exit()

#import astropy.io.fits as fits
#import os, sys, glob
#from os.path import join
##import pymangle as mangle
#import numpy as np
#import matplotlib.pyplot as p
#from scipy.interpolate import interp1d

#from astropy_healpix import healpy 

#import astropy.units as u
#import astropy.cosmology as cc
#cosmo = cc.Planck13
#from astropy.coordinates import SkyCoord

#from sklearn.neighbors import BallTree, DistanceMetric
#from astropy.table import Table,unique,Column
#from math import radians, cos, sin, asin, sqrt, pi

#import matplotlib
#matplotlib.use('Agg')
#matplotlib.rcParams.update({'font.size': 14})
#import matplotlib.pyplot as p

#agn_clustering_dir = '/data_XMM.data36s/comparat/AGN_clustering'
#agn_clustering_dir = '/home/comparat/data_XMM.data/AGN_clustering'

#catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )
#figure_dir = os.path.join(catalog_dir, 'figures_VAC' )

#if os.path.isdir(figure_dir)==False:
	#os.system('mkdir -p '+figure_dir)

#path_2_xmmsl_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_XMMSL2_DR16.fits')

#path_2_xmmsl_cat = os.path.join( catalog_dir, 'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits' )
#path_2_xmmsl_cat = os.path.join( catalog_dir, 'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits' )

#hd_xmmsl = fits.open(path_2_xmmsl_cat)[1].data_XMM.data

#prefix='XMMSL2_'

#data_XMM.data = hd_xmmsl

#ra = data_XMM.data['ALLW_RA']
#dec = data_XMM.data['ALLW_DEC']

#goodZ = (data_XMM.data['Z_BEST']>-1)

#print( len(data_XMM.data) )
#print( len(data_XMM.data[data_XMM.goodZ]), len(data_XMM.data[data_XMM.goodZ])*1./len(data_XMM.data) )
## fraction of target observed and fractino of successful redshifts as a function of quantities. 

#qty = 'XMMSL2_FLUX_B6'
#ok = data_XMM.data[qty]>0
#YY = np.log10(data_XMM.data[qty][ok])

#bins = np.arange(np.min(YY), np.max(YY), (np.max(YY)-np.min(YY))/20.)
#x_bins = (bins[1:]+bins[:-1])/2.

#out_all = np.histogram( YY, bins=bins )[0]
#out_zzz = np.histogram( YY[goodZ[ok]], bins=bins )[0]

#fig_out = os.path.join(figure_dir, prefix+qty+'_hist.png')

#p.figure(1, (5.5,5.5))
#p.axes([0.17, 0.17, 0.78, 0.75])
#p.tight_layout()
#p.plot(x_bins, out_all/np.sum(out_all), label='pdf')
##
#mid_line = out_zzz*1./out_all
#err_line = out_zzz**(-0.5)
#p.fill_between(x_bins, y1=mid_line*(1-err_line), y2=mid_line*(1+err_line), label='good Z/observed', alpha=0.5)
##
#p.xlabel(qty)#'Observed fraction per '+pixel_area_str+r' deg$^2$ pixels')
#p.ylabel('fraction')
#p.grid()
#p.ylim((0,1))
##p.yscale('log')
#p.legend(frameon=False, loc=4)
#p.savefig(fig_out)
#p.clf()

## NSIDE=2^5=32
## 
#nside_int = int(sys.argv[1])
#nside = 2**nside_int
#nside_str = 'hp'+str(nside_int)
#pixel_area = healpy.nside2pixarea(nside, degrees=True)
#pixel_area_str = str(np.round(pixel_area,2))

#hpx = healpy.ang2pix(nside, ra, dec, nest=True, lonlat=True)
#uniq_all = np.unique(hpx, return_counts=True)
#uniq_goodZ = np.unique(hpx[data_XMM.goodZ], return_counts=True)

#N_pixels_with_targets = len(uniq_all[0])
#print('N_pixels_with_targets', N_pixels_with_targets)
#N_pixels_with_goodZ = len(uniq_goodZ[0])
#print('N_pixels_with_goodZ', N_pixels_with_goodZ)

#all_targets = interp1d(uniq_all[0], uniq_all[1]*1.)
#success_rate = np.round( uniq_goodZ[1]/all_targets(uniq_goodZ[0]), 3)

#fig_out = os.path.join(figure_dir, prefix+'successZ_hist_'+nside_str+'.png')

#p.figure(1, (5.5,5.5))
#p.axes([0.17, 0.17, 0.78, 0.75])
#p.tight_layout()
##p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2)
#out = p.hist(success_rate, bins = np.arange(0.,0.99, 0.01, ), histtype='step', rasterized=True, lw=3)
##p.title(baseName)
#p.xlabel('Observed fraction per '+pixel_area_str+r' deg$^2$ pixels')
#p.ylabel('Number of pixels')
#p.grid()
##p.yscale('log')
##p.legend(frameon=False, loc=0)
#p.savefig(fig_out)
#p.clf()

#not_observed = np.in1d(uniq_all[0], uniq_goodZ[0], invert=True)

#ra_all, dec_all = healpy.pix2ang(nside, uniq_all[0][not_observed], nest=True, lonlat=True)
#ra, dec = healpy.pix2ang(nside, uniq_goodZ[0], nest=True, lonlat=True)

#fig_out = os.path.join(figure_dir, prefix+'successZ_ra_dec_'+nside_str+'.png')

#p.figure(1, (6.5,5.5))
#p.axes([0.17, 0.17, 0.78, 0.75])
#p.tight_layout()
##p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2)
#p.scatter(ra, dec, c=success_rate, s=14-2*nside_int, rasterized=True, cmap='Paired')
#p.colorbar(shrink=0.9)
#p.xlabel('R.A. [degree]')
#p.ylabel('Dec. [degree]')
#p.grid()
#p.title('Success rate per '+pixel_area_str+r' deg$^2$ pixels')
##p.yscale('log')
##p.legend(frameon=False, loc=0)
#p.savefig(fig_out)
#p.clf()

