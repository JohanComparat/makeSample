"""

python visual-inspection-stats.py

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

# on ds machines :
#agn_clustering_dir = '/data36s/comparat/AGN_clustering'
# on laptop
agn_clustering_dir = '/home/comparat/data/AGN_clustering'

catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )

git_dir = os.environ['GIT_MAKESAMPLE']
target_dir  = os.path.join(agn_clustering_dir, 'targets'  )
footprint_dir  = os.path.join(agn_clustering_dir, 'footprint'  )
figure_dir = os.path.join(git_dir, 'figures', 'agn', 'figures_VAC' )

if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)


path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_2RXS_DR16.fits')
VAC_2RXS = fits.open(path_2_cat)[1].data

path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_XMMSL2_DR16.fits')
VAC_XMMSL = fits.open(path_2_cat)[1].data

path_2_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits')
# contains 135,259
#hd0=fits.open(os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits'))
# contains 21,288
path_2_cat_dr16 = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')

# statistics of inspection :
#path_2_cat = os.path.join(catalog_dir, 'results_all.fits.gz')
#path_2_cat = os.path.join(catalog_dir, 'results_VAC_2RXS.fits')
data = VAC_2RXS # fits.open(path_2_cat)[1].data

s1 = (data['CONF_BEST']>=0)
def get_stats(s1):
	unique_id = (data['PLATE_BEST'][s1]*1e11 + data['MJD_BEST'][s1]*1e5 + data['FIBERID_BEST'][s1]).astype('int')
	uuu = np.unique(unique_id, return_index=True, return_inverse=True, return_counts=True)
	vis_stat = np.unique(uuu[3], return_counts=True)
	print(np.transpose([ vis_stat[0], vis_stat[1] ]))
	return uuu, vis_stat

print('0')
s1 = (data['CONF_BEST']==0)&(data['DR16_ZWARNING']>0)
get_stats(s1)

print('1')
s1 = (data['CONF_BEST']==1)&(data['DR16_ZWARNING']>0)
get_stats(s1)

print('2')
s1 = (data['CONF_BEST']==2)&(data['DR16_ZWARNING']>0)
get_stats(s1)

print('3')
s1 = (data['CONF_BEST']==3)&(data['DR16_ZWARNING']>0)
get_stats(s1)

# statistics of redshift measurement 



print('0')
s1 = (data['CONF_BEST']==0)&(data['DR16_ZWARNING']==0)
o0=get_stats(s1)

print('1')
s1 = (data['CONF_BEST']==1)&(data['DR16_ZWARNING']==0)
o1=get_stats(s1)

print('2')
s1 = (data['CONF_BEST']==2)&(data['DR16_ZWARNING']==0)
o2=get_stats(s1)

print('3')
s1 = (data['CONF_BEST']==3)&(data['DR16_ZWARNING']==0)
o3=get_stats(s1)


#print('0')
#s1 = (data['CONF_BEST']==0)&(data['DR16_ZWARNING']>0)
#o0=get_stats(s1)

#print('1')
#s1 = (data['CONF_BEST']==1)&(data['DR16_ZWARNING']>0)
#o1=get_stats(s1)

#print('2')
#s1 = (data['CONF_BEST']==2)&(data['DR16_ZWARNING']>0)
#o2=get_stats(s1)

#print('3')
#s1 = (data['CONF_BEST']==3)&(data['DR16_ZWARNING']>0)
#o3=get_stats(s1)




path_2_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits')
# contains 135,259
#hd0=fits.open(os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits'))
# contains 21,288
path_2_cat_dr16 = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')

MIN_SDSS_FIBER2MAG_i = 17.0
MAX_SDSS_FIBER2MAG_i = 22.5
MIN_SDSS_MODELMAG_i  = 16.0

EXI_ML_min_str = sys.argv[1]

if EXI_ML_min_str=='6p5':
	EXI_ML_min=6.5
if EXI_ML_min_str=='10':
	EXI_ML_min=10

prefix='2RXS_ExiML'+EXI_ML_min_str+'_'

#prefix='XMMSL2_ExiML10_'
#prefix='XMMSL2_ExiML6p5_'

#def get_arrays(path_2_cat, EXI_ML_min = 10.):
def get_arrays_obs(path_2_cat, EXI_ML_min = EXI_ML_min):
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
	c1 = (observed) & ((data['Z_BEST']>0) | ((data['DR16_Z']>0) & (data['DR16_Z_ERR']>0)))
	c2 = (c1) & (data['CONF_BEST']==3)
	c3 = (c1) & (data['CONF_BEST']==2) & ((data['CLASS_BEST']=='BLAZAR')|(data['CLASS_BEST']=='BLLAC'))
	c4 = (c1) & (data['DR16_SN_MEDIAN_ALL']>=2) & (data['DR16_ZWARNING']==0 )
	c5 = (c1) & (data['VI_REINSPECT_FLAG'] == 0) & (data['VI_NINSPECTORS']>2)
	c6 = (c1) & (data['VI_AM_RECONCILED']==1)
	goodZ =  (c2) | (c3) | (c4) | (c5) | (c6) 
	#goodZ     = (observed) & ( 
		#(data['CONF_BEST']==3) | 
		#( (data['CONF_BEST']>=2) & 
			#( (data['CLASS_BEST']=='BLAZAR') | (data['CLASS_BEST']=='BLLAC') )  
		#) )
	print(len(ra), len(ra[(high_conf)]))
	return data[observed], goodZ[observed]

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
	c1 = (observed) & ((data['Z_BEST']>0) | ((data['DR16_Z']>0) & (data['DR16_Z_ERR']>0)))
	c2 = (c1) & (data['CONF_BEST']==3)
	c3 = (c1) & (data['CONF_BEST']==2) & ((data['CLASS_BEST']=='BLAZAR')|(data['CLASS_BEST']=='BLLAC'))
	c4 = (c1) & (data['DR16_SN_MEDIAN_ALL']>=2) & (data['DR16_ZWARNING']==0 )
	c5 = (c1) & (data['VI_REINSPECT_FLAG'] == 0) & (data['VI_NINSPECTORS']>2)
	c6 = (c1) & (data['VI_AM_RECONCILED']==1)
	goodZ =  (c2) | (c3) | (c4) | (c5) | (c6) 
	#goodZ     = (observed) & ( 
		#(data['CONF_BEST']==3) | 
		#( (data['CONF_BEST']>=2) & 
			#( (data['CLASS_BEST']=='BLAZAR') | (data['CLASS_BEST']=='BLLAC') )  
		#) )
	print(len(ra), len(ra[high_conf]) )
	return data[(high_conf)], ra[(high_conf)], dec[(high_conf)], targeted[(high_conf)], observed[(high_conf)], goodZ[(high_conf)], hd_rass_i


data, goodZ = get_arrays_obs(path_2_cat_dr16)

print(data['DR16_Z'], data['DR16Q_Z'], data['Z_BEST'])


cb0=(data['CONF_BEST']==0)
cb1=(data['CONF_BEST']==1)
cb2=(data['CONF_BEST']==2)
cb3=(data['CONF_BEST']==3)
cb9=(data['CONF_BEST']==-99)

print('0', np.unique(data['CLASS_BEST'][cb0], return_counts=True))
print('1', np.unique(data['CLASS_BEST'][cb1], return_counts=True))
print('2', np.unique(data['CLASS_BEST'][cb2], return_counts=True))
print('3', np.unique(data['CLASS_BEST'][cb3], return_counts=True))
print('-99', np.unique(data['CLASS_BEST'][cb9], return_counts=True))
out_name, out_N = np.unique(data['CLASS_BEST'][goodZ], return_counts=True)
for aa,bb in zip(out_name, out_N):
	print(aa, bb )

qty='SDSS_FIBER2MAG_i'

good = (data[qty]>0)
goodZED = (good)&(goodZ)
cb0=(good)&(data['CONF_BEST']==0)
cb1=(good)&(data['CONF_BEST']==1)
cb2=(good)&(data['CONF_BEST']==2)
cb3=(good)&(data['CONF_BEST']==3)
cb9=(good)&(data['CONF_BEST']==-99)

fig_out = os.path.join(figure_dir, prefix+'_'+qty+'_hist_dr16.png')

YY = data[qty][good]#'DR16Q_SN_MEDIAN_ALL']

XMIN = 16 
XMAX = 30
DX= 0.5
bins = np.arange(XMIN, XMAX+DX, DX)# (np.max(YY)-np.min(YY))/10.)
x_bins = (bins[1:]+bins[:-1])/2.

p.figure(1, (5.5,4.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()

out_all = np.histogram( data[qty][good], bins=bins)[0]
out_cb0 = np.histogram( data[qty][ cb0], bins=bins)[0]
out_cb1 = np.histogram( data[qty][ cb1], bins=bins)[0]
out_cb2 = np.histogram( data[qty][ cb2], bins=bins)[0]
out_cb3 = np.histogram( data[qty][ cb3], bins=bins)[0]

#p.plot(x_bins, out_all/np.sum(out_all),lw=3, label='pdf')
#p.plot(x_bins, out_cb0/np.sum(out_cb0),lw=2, label='0')
p.hist(data[qty][good], bins=bins, histtype='step', lw=4, label='all')#, N='+str(int(np.sum(out_all))) )
#p.hist(data[qty][goodZED], bins=bins, histtype='step', lw=2, label='goodZ', ls='dashed')#, N='+str(int(np.sum(out_all))) )
p.hist(data[qty][ cb3], bins=bins, histtype='step',lw=2, label='3')#, N='+str(int(np.sum(out_cb3))), ls='dashed' )
p.hist(data[qty][ cb2], bins=bins, histtype='step',lw=3, label='2')#, N='+str(int(np.sum(out_cb2))), ls='dashed' )
p.hist(data[qty][ cb1], bins=bins, histtype='step',lw=3, label='1')#, N='+str(int(np.sum(out_cb1))) )
p.hist(data[qty][ cb0], bins=bins, histtype='step',lw=3, label='0')#, N='+str(int(np.sum(out_cb0))) )
p.hist(data[qty][ cb9], bins=bins, histtype='step',lw=3, label='-99')#, N='+str(int(np.sum(out_cb0))) )
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
p.ylabel('N')
p.grid()
p.xlim((16.8, 22.7 ))# np.max(x_bins[out_all>1])  ))
#p.ylim((0.9,5000))
p.yscale('log')
#p.legend(frameon=True, loc=6)
p.savefig(fig_out)
p.clf()

qty='SDSS_MODELMAG_i'

good = (data[qty]>0)
goodZED = (good)&(goodZ)
cb0=(good)&(data['CONF_BEST']==0)
cb1=(good)&(data['CONF_BEST']==1)
cb2=(good)&(data['CONF_BEST']==2)
cb3=(good)&(data['CONF_BEST']==3)
cb9=(good)&(data['CONF_BEST']==-99)

fig_out = os.path.join(figure_dir, prefix+'_'+qty+'_hist_dr16.png')

YY = data[qty][good]#'DR16Q_SN_MEDIAN_ALL']

XMIN = 16 
XMAX = 30
DX= 0.5 # 1.0 #0.5
bins = np.arange(XMIN, XMAX+DX, DX)# (np.max(YY)-np.min(YY))/10.)
x_bins = (bins[1:]+bins[:-1])/2.

p.figure(1, (5.5,4.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()

out_all = np.histogram( data[qty][good], bins=bins)[0]
out_cb0 = np.histogram( data[qty][ cb0], bins=bins)[0]
out_cb1 = np.histogram( data[qty][ cb1], bins=bins)[0]
out_cb2 = np.histogram( data[qty][ cb2], bins=bins)[0]
out_cb3 = np.histogram( data[qty][ cb3], bins=bins)[0]
out_cb9 = np.histogram( data[qty][ cb9], bins=bins)[0]

#p.plot(x_bins, out_all/np.sum(out_all),lw=3, label='pdf')
p.hist(data[qty][good], bins=bins, histtype='step', lw=4, label='all')#, N='+str(int(np.sum(out_all))) )
#p.hist(data[qty][goodZED], bins=bins, histtype='step', lw=2, label='goodZ', ls='dashed')#, N='+str(int(np.sum(out_all))) )
p.hist(data[qty][ cb3], bins=bins, histtype='step',lw=2, label='3')#, N='+str(int(np.sum(out_cb3))), ls='dashed' )
p.hist(data[qty][ cb2], bins=bins, histtype='step',lw=3, label='2')#, N='+str(int(np.sum(out_cb2))), ls='dashed' )
p.hist(data[qty][ cb1], bins=bins, histtype='step',lw=3, label='1')#, N='+str(int(np.sum(out_cb1))) )
p.hist(data[qty][ cb0], bins=bins, histtype='step',lw=3, label='0')#, N='+str(int(np.sum(out_cb0))) )
p.hist(data[qty][ cb9], bins=bins, histtype='step',lw=3, label='-99')#, N='+str(int(np.sum(out_cb0))) )
#p.plot(x_bins, out_cb0/np.sum(out_cb0),lw=3, label='0, N='+str(int(np.sum(out_cb0))) )
#p.plot(x_bins, out_cb1/np.sum(out_cb1),lw=3, label='1, N='+str(int(np.sum(out_cb1))) )
#p.plot(x_bins, out_cb2/np.sum(out_cb2),lw=3, label='2, N='+str(int(np.sum(out_cb2))), ls='dashed'  )
#p.plot(x_bins, out_cb3/np.sum(out_cb3),lw=3, label='3, N='+str(int(np.sum(out_cb3))), ls='dashed'  )
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
p.ylabel('N')
p.grid()
p.xlim((15.8, 22.2 )) #np.max(x_bins[out_all>1]) + DX ))
#p.ylim((0,1))
p.yscale('log')
p.legend(frameon=True, loc=2)
p.savefig(fig_out)
p.clf()


qty='DR16Q_SN_MEDIAN_ALL'

good = (data[qty]>0)
goodZED = (good)&(goodZ)
cb0=(good)&(data['CONF_BEST']==0)
cb1=(good)&(data['CONF_BEST']==1)
cb2=(good)&(data['CONF_BEST']==2)
cb3=(good)&(data['CONF_BEST']==3)
cb9=(good)&(data['CONF_BEST']==-99)

fig_out = os.path.join(figure_dir, prefix+'_'+qty+'_hist_dr16.png')

YY = data[qty][good]#'DR16Q_SN_MEDIAN_ALL']

XMIN = -1 
XMAX = 2
DX=0.25
bins = 10**np.arange(XMIN, XMAX+DX, DX)# (np.max(YY)-np.min(YY))/10.)
x_bins = (bins[1:]+bins[:-1])/2.

p.figure(1, (5.5,4.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()

out_all = np.histogram( data[qty][good], bins=bins)[0]
out_cb0 = np.histogram( data[qty][ cb0], bins=bins)[0]
out_cb1 = np.histogram( data[qty][ cb1], bins=bins)[0]
out_cb2 = np.histogram( data[qty][ cb2], bins=bins)[0]
out_cb3 = np.histogram( data[qty][ cb3], bins=bins)[0]

p.hist(data[qty][good], bins=bins, histtype='step', lw=4, label='all')#, N='+str(int(np.sum(out_all))) )
#p.hist(data[qty][goodZED], bins=bins, histtype='step', lw=2, label='goodZ', ls='dashed')#, N='+str(int(np.sum(out_all))) )
p.hist(data[qty][ cb3], bins=bins, histtype='step',lw=2, label='3')#, N='+str(int(np.sum(out_cb3))), ls='dashed' )
p.hist(data[qty][ cb2], bins=bins, histtype='step',lw=3, label='2')#, N='+str(int(np.sum(out_cb2))), ls='dashed' )
p.hist(data[qty][ cb1], bins=bins, histtype='step',lw=3, label='1')#, N='+str(int(np.sum(out_cb1))) )
p.hist(data[qty][ cb0], bins=bins, histtype='step',lw=3, label='0')#, N='+str(int(np.sum(out_cb0))) )
p.hist(data[qty][ cb9], bins=bins, histtype='step',lw=3, label='-99')#, N='+str(int(np.sum(out_cb0))) )

#p.plot(x_bins, out_all/np.sum(out_all),lw=3, label='pdf')
##p.plot(x_bins, out_cb0/np.sum(out_cb0),lw=2, label='0')
#p.plot(x_bins, out_cb1/np.sum(out_cb1),lw=2, label='1, N='+str(int(np.sum(out_cb1))) )
#p.plot(x_bins, out_cb2/np.sum(out_cb2),lw=2, label='2, N='+str(int(np.sum(out_cb2))) )
#p.plot(x_bins, out_cb3/np.sum(out_cb3),lw=2, label='3, N='+str(int(np.sum(out_cb3))) )
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

p.xlabel(r'DR16Q_SN_MEDIAN_ALL')
p.ylabel('N')
p.grid()
#p.xlim((bins.min(), bins.max()))
p.xlim((0.4, 150))
#p.ylim((0,1))
p.xscale('log')
p.yscale('log')
#p.legend(frameon=True, loc=0)
p.savefig(fig_out)
p.clf()
