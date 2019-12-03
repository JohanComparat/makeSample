"""

python dr16-completeness-VAC-2RXS.py 4 6p5
python dr16-completeness-VAC-2RXS.py 4 10

"""
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

EXI_ML_min_str = sys.argv[2]

if EXI_ML_min_str=='6p5':
	EXI_ML_min=6.5
if EXI_ML_min_str=='10':
	EXI_ML_min=10

prefix='2RXS_ExiML'+EXI_ML_min_str+'_'
# prefix='XMMSL2_ExiML10_'
# prefix='XMMSL2_ExiML6p5_'
# 2RXS
qty = 'RXS_SRC_FLUX'
# XMMSL2
# qty = 'XMMSL2_FLUX_B6'

class Catalog: pass

data, high_conf, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, bl, nl, class_best = get_arrays_2rxs(path_2_cat_2RXS = path_2_cat_2RXS, EXI_ML_min = EXI_ML_min)

data_RXS = Catalog()
data_RXS.data = data
data_RXS.high_conf = high_conf
data_RXS.targeted = targeted
data_RXS.observed = observed
data_RXS.goodZ = goodZ
data_RXS.idZ = idZ
data_RXS.agnZ = agnZ
data_RXS.clusters = clusters
data_RXS.blazars_noZ = blazars_noZ
data_RXS.stars = stars
data_RXS.bl = bl
data_RXS.nl = nl
data_RXS.class_best = class_best

fig_out = os.path.join(figure_dir, prefix+qty+'_hist_dr16.png')

XMIN = -13. 
XMAX = -10
DX = 1.
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
#print(len(all_data_RXS.data[(all_data_RXS.data['RXS_IN_BOSS']==1)])) 
print('All         &', len(data_RXS.data) )
print('All         &', len(data_RXS.data[data_RXS.high_conf]) )
print('Targets     &', len(data_RXS.data[data_RXS.targeted]) )#, len(data_RXS.data[data_RXS.targeted])*1./len(data_RXS.data) )
print('Observed    &', len(data_RXS.data[data_RXS.observed]) )#, len(data_RXS.data[data_RXS.observed])*1./len(data_RXS.data[data_RXS.targeted]) )
print('Identified  &', len(data_RXS.data[data_RXS.idZ])      )#, len(data_RXS.data[data_RXS.idZ])*1./len(data_RXS.data[data_RXS.observed]) )
print('Blazar no Z &', len(data_RXS.data[data_RXS.blazars_noZ]))#, len(data_RXS.data[data_RXS.blazars_noZ])*1./len(data_RXS.data[data_RXS.observed]) )
print('good Z      &', len(data_RXS.data[data_RXS.goodZ])    )#, len(data_RXS.data[data_RXS.goodZ])*1./len(data_RXS.data[data_RXS.observed]) )
print('AGN         &', len(data_RXS.data[data_RXS.agnZ])     )#, len(data_RXS.data[data_RXS.agnZ])*1./len(data_RXS.data[data_RXS.observed]) )
print('Cluster     &', len(data_RXS.data[data_RXS.clusters]) )#, len(data_RXS.data[data_RXS.clusters])*1./len(data_RXS.data[data_RXS.observed]) )
print('Star        &', len(data_RXS.data[data_RXS.stars])    )#, len(data_RXS.data[data_RXS.stars])*1./len(data_RXS.data[data_RXS.observed]) )
# fraction of target observed and fractino of successful redshifts as a function of quantities. 
print('===============')
print( len(data_RXS.data[data_RXS.observed]), len(data_RXS.data[data_RXS.observed])*1./len(data_RXS.data[data_RXS.observed]) )
out = np.unique(data_RXS.class_best[data_RXS.clusters], return_counts=True)
print(out)

#cb = data_RXS.data['CONF_BEST']
#for cb_val in np.unique(cb):
	#cb_s = (cb==cb_val)
	#print(cb_val, ' & ', len(data_RXS.data[observed & cb_s]) )

YY = np.log10(data_RXS.data[qty])

out_all = np.histogram( YY, bins=bins )[0]
out_tar = np.histogram( YY[data_RXS.targeted], bins=bins )[0]
out_obs = np.histogram( YY[data_RXS.observed], bins=bins )[0]
out_zzz = np.histogram( YY[data_RXS.goodZ], bins=bins )[0]

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

#p.plot(x_bins, out_all/np.sum(out_all), label='pdf')
p.hist(YY, bins=bins, label='normed hist', weights = np.ones_like(YY)/np.sum(out_all), histtype='step', lw=3)

p.xlabel(r'$\log_{10}(F_X [erg.\, cm^{-2}. s^{-1}])$')
p.ylabel('fraction')
p.grid()
p.xlim((XMIN, XMAX))
p.ylim((-0.02,1.02))
#p.yscale('log')
p.legend(frameon=True, loc=3)
p.savefig(fig_out)
p.clf()

sys.exit()

#ra  = data_RXS.data['ALLW_RAJ2000']
#dec = data_RXS.data['ALLW_DEJ2000']

#hpx = healpy.ang2pix(nside, ra, dec, nest=True, lonlat=True)
#uniq_all = np.unique(hpx, return_counts=True)
#uniq_targeted = np.unique(hpx[data_RXS.targeted], return_counts=True)
#uniq_observed = np.unique(hpx[data_RXS.observed], return_counts=True)
#uniq_goodZ = np.unique(hpx[data_RXS.goodZ], return_counts=True)

#N_pixels_with_targets = len(uniq_all[0])
#print('N_pixels_with_targets', N_pixels_with_targets)
#N_pixels_with_observed_targets = len(uniq_observed[0])
#print('N_pixels_with_observed_targets', N_pixels_with_observed_targets)
#N_pixels_with_goodZ = len(uniq_goodZ[0])
#print('N_pixels_with_goodZ', N_pixels_with_goodZ)

#all_targets = interp1d(uniq_targeted[0], uniq_targeted[1]*1.)
#completeness = np.round( uniq_observed[1]/all_targets(uniq_observed[0]), 3)
#success_rate = np.round( uniq_goodZ[1]/all_targets(uniq_goodZ[0]), 3)

#fig_out = os.path.join(figure_dir, prefix+'successZ_hist_'+nside_str+'.png')

#p.figure(1, (5.5,5.5))
#p.axes([0.17, 0.17, 0.78, 0.75])
#p.tight_layout()
#out1 = p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2, label='O/T')
#out  = p.hist(success_rate, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=3, label='Z/O')
##p.title(baseName)
#p.xlabel('Observed fraction per '+pixel_area_str+r' deg$^2$ pixels')
#p.ylabel('Number of pixels')
#p.grid()
##p.yscale('log')
#p.legend(frameon=False, loc=0)
#p.savefig(fig_out)
#p.clf()

#not_observed = np.in1d(uniq_all[0], uniq_observed[0], invert=True)

#ra_all, dec_all = healpy.pix2ang(nside, uniq_all[0][not_observed], nest=True, lonlat=True)
#ra_hp, dec_hp = healpy.pix2ang(nside, uniq_observed[0], nest=True, lonlat=True)

#p2map = os.path.join(catalog_dir, 'DR16_footprint_hpx10.fits')
#p2map = os.path.join(catalog_dir, 'DR16_wSEQUELS_footprint_clipped_complete_hpx10.fits')
#p2map = os.path.join(catalog_dir, 'DR16_wSEQUELS_footprint_hpx10.fits')
#hd=fits.open(p2map)

#NSIDE = 1024
#all_pix = np.arange(healpy.nside2npix(NSIDE))
#ra_mask, dec_mask = healpy.pix2ang(NSIDE, all_pix[hd[1].data_RXS.data['MASK'].astype('bool')], nest=True, lonlat=True)

#fig_out = os.path.join(figure_dir, prefix+'successZ_ra_dec_'+nside_str+'.png')

#p.figure(1, (7.,5.5))
#p.axes([0.17, 0.17, 0.78, 0.75])
#p.tight_layout()
##p.plot(ra_mask, dec_mask, 'ko', rasterized=True, alpha=0.0015)
##p.hist(completeness, bins = np.arange(0.,1.1, 0.05, ), histtype='step', rasterized=True, lw=2)
#p.scatter(ra_hp, dec_hp, c=completeness, s=20, marker='o', rasterized=True, cmap='RdBu', label='exiML>'+str(np.round(EXI_ML_min,1)))#, alpha=0.75)
#p.colorbar(shrink=0.9)
#p.xlabel('R.A. [degree]')
#p.ylabel('Dec. [degree]')
#p.grid()
#p.title('Fraction observed per '+pixel_area_str+r' deg$^2$ pixels')
##p.yscale('log')
#p.legend(frameon=False, loc=0)
#p.savefig(fig_out)
#p.clf()

#fig_out = os.path.join(figure_dir, 'eBOSS-dr16-wSequels.png')

#p.figure(1, (6.5,5.5))
#p.axes([0.17, 0.17, 0.78, 0.75])
#p.tight_layout()
#p.plot(ra_mask, dec_mask, 'ko', rasterized=True)
#p.xlabel('R.A. [degree]')
#p.ylabel('Dec. [degree]')
#p.grid()
##p.yscale('log')
##p.legend(frameon=False, loc=0)
#p.savefig(fig_out)
#p.clf()


