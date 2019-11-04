
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

# on ds machines :
#agn_clustering_dir = '/data36s/comparat/AGN_clustering'
# on laptop
agn_clustering_dir = '/home/comparat/data/AGN_clustering'

nside_int = int(sys.argv[1])
nside = 2**nside_int
nside_str = 'hp'+str(nside_int)
pixel_area = healpy.nside2pixarea(nside, degrees=True)
pixel_area_str = str(np.round(pixel_area,2))

catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )

git_dir = os.environ['GIT_MAKESAMPLE']
target_dir  = os.path.join(agn_clustering_dir, 'targets'  )
footprint_dir  = os.path.join(agn_clustering_dir, 'footprint'  )
figure_dir = os.path.join(git_dir, 'figures', 'agn', 'figures_VAC_XMMSL' )

if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)

#path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_2RXS_DR16.fits')

#path_2_2RXS_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits') # '2RXS_AllWISE_catalog_paper_2017May26_v5_11_0_sdss_26_VERON_MASKED_GAIA_star_mask.fits')
#full_2RXS = fits.open(path_2_2RXS_cat)[1].data

#path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_2RXS_DR16.fits')
#VAC_2RXS = fits.open(path_2_cat)[1].data

path_2_cat_xmmsl = os.path.join(catalog_dir, 'VAC_SPIDERS_XMMSL2_DR16.fits')
VAC_XMMSL = fits.open(path_2_cat_xmmsl)[1].data

path_2_cat_xmmsl = os.path.join(catalog_dir, 'XMMSL2',  'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_REDMAPPER.fits')
VAC_XMMSL = fits.open(path_2_cat_xmmsl)[1].data

#path_2_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits')
# contains 135,259
#hd0=fits.open(os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits'))
# contains 21,288
#path_2_cat_dr16 = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')
#path_2_cat_dr16 = os.path.join(catalog_dir, '2RXS', 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER.fits')

MIN_SDSS_FIBER2MAG_i = 17.0
MAX_SDSS_FIBER2MAG_i = 22.5
MIN_SDSS_MODELMAG_i  = 16.0

EXI_ML_min_str = '10'
EXI_ML_min = 10

prefix='XMMSL2_ExiML'+EXI_ML_min_str+'_'

#prefix='XMMSL2_ExiML10_'
#prefix='XMMSL2_ExiML6p5_'

#def get_arrays(path_2_cat, EXI_ML_min = 10.):
def get_arrays_xmmsl2(path_2_cat, EXI_ML_min = EXI_ML_min):
	hd_rass_i = fits.open(path_2_cat)[1].data
	targets = (hd_rass_i['NWAY_match_flag']!=2) & (hd_rass_i['XMMSL2_IN_BOSS']==1) & (hd_rass_i['FLAG_SDSSv5b_best']==1) & (hd_rass_i['NWAY_p_any']>=0.01) & (hd_rass_i['NUM_SDSS']>=1)
	#  NWAY_match_flag!=2 &&RXS_IN_BOSS==1 &&FLAG_SDSSv5b_best==1 &&
	hd_rass = fits.open(path_2_cat)[1].data[targets]
	#print(np.unique(hd_rass['CLASS_BEST']))
	data = hd_rass
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
	print(len(ra), len(ra[high_conf]) )
	return data[(high_conf)], ra[(high_conf)], dec[(high_conf)], targeted[(high_conf)], observed[(high_conf)], goodZ[(high_conf)], idZ[(high_conf)], agnZ[(high_conf)], clusters[(high_conf)], blazars_noZ[(high_conf)], stars[(high_conf)], hd_rass_i

data, ra, dec, targeted, observed, goodZ, idZ, agnZ, clusterZ, blazarNOZ, stars, all_data = get_arrays_xmmsl2(path_2_cat_xmmsl)

#data['RXS_ExtML'    ]
#data['RXS_ExiML1RXS'] 
#data['RXS_ExiML2RXS']

# 2RXS
#qty = 'RXS_SRC_FLUX'
# XMMSL2
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
	data['XMMSL2_DET_ML_B6'] , 
	data['XMMSL2_DET_ML_B7'] , 
	data['XMMSL2_DET_ML_B8'] ])
ExiML[np.isnan(ExiML)]=0
ExiML_arr = np.max(ExiML,axis=1) 


print(len(all_data[(all_data['XMMSL2_IN_BOSS']==1)])) 
print('All         &', len(data) )
print('Targets     &', len(data[targeted]) )#, len(data[targeted])*1./len(data) )
print('Observed    &', len(data[observed]) )#, len(data[observed])*1./len(data[targeted]) )
print('Identified  &', len(data[idZ])      )#, len(data[idZ])*1./len(data[observed]) )
print('Blazar no Z &', len(data[blazarNOZ]))#, len(data[blazarNOZ])*1./len(data[observed]) )
print('good Z      &', len(data[goodZ])    )#, len(data[goodZ])*1./len(data[observed]) )
print('AGN         &', len(data[agnZ])     )#, len(data[agnZ])*1./len(data[observed]) )
print('Cluster     &', len(data[clusterZ]) )#, len(data[clusterZ])*1./len(data[observed]) )
print('Star        &', len(data[stars])    )#, len(data[stars])*1./len(data[observed]) )

# fraction of target observed and fractino of successful redshifts as a function of quantities. 

#print('===============')
#print( len(data[observed]), len(data[observed])*1./len(data[observed]) )

#cb = data['CONF_BEST']
#for cb_val in np.unique(cb):
	#cb_s = (cb==cb_val)
	#print(cb_val, ' & ', len(data[observed & cb_s]) )#, len(data[observed & cb_s])*1./len(data[observed]) )


#YY = np.log10(data[qty])

#out_all = np.histogram( YY, bins=bins )[0]
#out_tar = np.histogram( YY[targeted], bins=bins )[0]
#out_obs = np.histogram( YY[observed], bins=bins )[0]
#out_zzz = np.histogram( YY[goodZ], bins=bins )[0]

##
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

##p.plot(x_bins, out_all/np.sum(out_all), label='pdf')
#p.hist(YY, bins=bins, label='normed hist', weights = np.ones_like(YY)/np.sum(out_all), histtype='step', lw=3)

#p.xlabel(r'$\log_{10}(F_X [erg.\, cm^{-2}. s^{-1}])$')
#p.ylabel('fraction')
#p.grid()
#p.xlim((XMIN, XMAX))
#p.ylim((-0.02,1.02))
##p.yscale('log')
#p.legend(frameon=True, loc=3)
#p.savefig(fig_out)
#p.clf()


dz = 0.1
z_min=0.0005 
redshift = data['DR16_Z']
z_max = np.max(redshift[agnZ]) # np.max([np.max(redshift),np.max(Z_XMMSL2)])
zs = np.arange(z_min, z_max + dz, dz)
xzs = (zs[1:]+zs[:-1])/2.

NN_xmmsl2 = np.histogram(redshift[agnZ], bins=zs)[0]
#density = NN_xmmsl2/area_eboss_dr16#/np.max(NN)

NN_clusters = np.histogram(redshift[clusterZ], bins=zs)[0]
#density = NN_clusters/area_eboss_dr16#/np.max(NN)

np.savetxt(os.path.join( figure_dir, "histogram_redshift.txt"), np.transpose([zs[:-1], zs[1:], NN_xmmsl2]), delimiter=" & ", newline=" \\\\ \n", fmt='%10.1f')

out = np.unique(data['CLASS_BEST'][clusterZ], return_counts=True)
print(out)

agn_in_cluster = (clusterZ) & (data['CLASS_BEST']!="GALAXY")
for el in np.transpose([data['ALLW_RA'][agn_in_cluster],data['ALLW_DEC'][agn_in_cluster], data['DR16_Z'][agn_in_cluster], data['REDMAPPER_Z_LAMBDA'][agn_in_cluster], data['REDMAPPER_Separation'][agn_in_cluster]]):
	print(np.round(el,5))


sys.exit()

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
figure_dir = os.path.join(catalog_dir, 'figures_VAC' )

if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)

path_2_xmmsl_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_XMMSL2_DR16.fits')

path_2_xmmsl_cat = os.path.join( catalog_dir, 'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits' )
path_2_xmmsl_cat = os.path.join( catalog_dir, 'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits' )

hd_xmmsl = fits.open(path_2_xmmsl_cat)[1].data

prefix='XMMSL2_'

data = hd_xmmsl

ra = data['ALLW_RA']
dec = data['ALLW_DEC']

goodZ = (data['Z_BEST']>-1)

print( len(data) )
print( len(data[goodZ]), len(data[goodZ])*1./len(data) )
# fraction of target observed and fractino of successful redshifts as a function of quantities. 

qty = 'XMMSL2_FLUX_B6'
ok = data[qty]>0
YY = np.log10(data[qty][ok])

bins = np.arange(np.min(YY), np.max(YY), (np.max(YY)-np.min(YY))/20.)
x_bins = (bins[1:]+bins[:-1])/2.

out_all = np.histogram( YY, bins=bins )[0]
out_zzz = np.histogram( YY[goodZ[ok]], bins=bins )[0]

fig_out = os.path.join(figure_dir, prefix+qty+'_hist.png')

p.figure(1, (5.5,5.5))
p.axes([0.17, 0.17, 0.78, 0.75])
p.tight_layout()
p.plot(x_bins, out_all/np.sum(out_all), label='pdf')
#
mid_line = out_zzz*1./out_all
err_line = out_zzz**(-0.5)
p.fill_between(x_bins, y1=mid_line*(1-err_line), y2=mid_line*(1+err_line), label='good Z/observed', alpha=0.5)
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

hpx = healpy.ang2pix(nside, ra, dec, nest=True, lonlat=True)
uniq_all = np.unique(hpx, return_counts=True)
uniq_goodZ = np.unique(hpx[goodZ], return_counts=True)

N_pixels_with_targets = len(uniq_all[0])
print('N_pixels_with_targets', N_pixels_with_targets)
N_pixels_with_goodZ = len(uniq_goodZ[0])
print('N_pixels_with_goodZ', N_pixels_with_goodZ)

all_targets = interp1d(uniq_all[0], uniq_all[1]*1.)
success_rate = np.round( uniq_goodZ[1]/all_targets(uniq_goodZ[0]), 3)

fig_out = os.path.join(figure_dir, prefix+'successZ_hist_'+nside_str+'.png')

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

not_observed = np.in1d(uniq_all[0], uniq_goodZ[0], invert=True)

ra_all, dec_all = healpy.pix2ang(nside, uniq_all[0][not_observed], nest=True, lonlat=True)
ra, dec = healpy.pix2ang(nside, uniq_goodZ[0], nest=True, lonlat=True)

fig_out = os.path.join(figure_dir, prefix+'successZ_ra_dec_'+nside_str+'.png')

p.figure(1, (6.5,5.5))
p.axes([0.17, 0.17, 0.78, 0.75])
p.tight_layout()
#p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2)
p.scatter(ra, dec, c=success_rate, s=14-2*nside_int, rasterized=True, cmap='Paired')
p.colorbar(shrink=0.9)
p.xlabel('R.A. [degree]')
p.ylabel('Dec. [degree]')
p.grid()
p.title('Success rate per '+pixel_area_str+r' deg$^2$ pixels')
#p.yscale('log')
#p.legend(frameon=False, loc=0)
p.savefig(fig_out)
p.clf()

