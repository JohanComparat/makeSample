"""

python visual-inspection-stats.py

"""
from catalog_lib import *

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

#data_RXS.class_best

#data_XMM.class_best
"""
identical = (data_XMM.data['CLASS_BEST'] == data_XMM.data['merged_class'])
DATA = np.transpose([ data_XMM.data['PLATE_BEST'] , data_XMM.data['MJD_BEST'] , data_XMM.data['FIBERID_BEST'],  data_XMM.data['merged_class'] ])[identical==False]
for el in DATA:
	print(" & ".join(el)+" \\\\")

identical = (data_RXS.data['CLASS_BEST'] == data_RXS.data['merged_class'])
DATA = np.transpose([ data_RXS.data['PLATE_BEST'] , data_RXS.data['MJD_BEST'] , data_RXS.data['FIBERID_BEST'],  data_RXS.data['merged_class'] ])[identical==False]
for el in DATA:
	print(" & ".join(el)+" \\\\")


identical = (data_XMM.data['CLASS_BEST'] == data_XMM.data['merged_class']) & (data_XMM.data['CLASS_BEST']=="STAR") & (data_XMM.data['XID_XID']>0)
DATA = np.transpose([ data_XMM.data['PLATE_BEST'] , data_XMM.data['MJD_BEST'] , data_XMM.data['FIBERID_BEST'],  data_XMM.data['XID_XID'],  data_XMM.data['XID_SUBCLASS'] ])[identical]
for el in DATA[::-1]:
	print(" & ".join(el)+" \\\\")

identical = (data_RXS.data['CLASS_BEST'] == data_RXS.data['merged_class'])& (data_RXS.data['CLASS_BEST']=="STAR") & (data_RXS.data['XID_XID']>0)
DATA = np.transpose([ data_RXS.data['PLATE_BEST'] , data_RXS.data['MJD_BEST'] , data_RXS.data['FIBERID_BEST'],  data_RXS.data['XID_XID'],  data_RXS.data['XID_SUBCLASS'] ])[identical]
ids = np.argsort(data_RXS.data['PLATE_BEST'][identical])
for el in DATA[ids]:
	print(" & ".join(el)+" \\\\")
"""
#path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_2RXS_DR16.fits')
#VAC_2RXS = fits.open(path_2_cat)[1].data

#path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_XMMSL2_DR16.fits')
#VAC_XMMSL = fits.open(path_2_cat)[1].data

#path_2_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits')
## contains 135,259
## hd0=fits.open(os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits'))
## contains 21,288
#path_2_cat_dr16 = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')

## statistics of inspection :
##path_2_cat = os.path.join(catalog_dir, 'results_all.fits.gz')
##path_2_cat = os.path.join(catalog_dir, 'results_VAC_2RXS.fits')
#data = VAC_2RXS # fits.open(path_2_cat)[1].data

#s1 = (data['CONF_BEST']>=0)
#def get_stats(s1):
	#unique_id = (data['PLATE_BEST'][s1]*1e11 + data['MJD_BEST'][s1]*1e5 + data['FIBERID_BEST'][s1]).astype('int')
	#uuu = np.unique(unique_id, return_index=True, return_inverse=True, return_counts=True)
	#vis_stat = np.unique(uuu[3], return_counts=True)
	#print(np.transpose([ vis_stat[0], vis_stat[1] ]))
	#return uuu, vis_stat

# statistics of redshift measurement 
#print('0')
#s1 = (data['CONF_BEST']==0)&(data['DR16_ZWARNING']>0)
#get_stats(s1)

#print('1')
#s1 = (data['CONF_BEST']==1)&(data['DR16_ZWARNING']>0)
#get_stats(s1)

#print('2')
#s1 = (data['CONF_BEST']==2)&(data['DR16_ZWARNING']>0)
#get_stats(s1)

#print('3')
#s1 = (data['CONF_BEST']==3)&(data['DR16_ZWARNING']>0)
#get_stats(s1)

#path_2_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits')
# contains 135,259
#hd0=fits.open(os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits'))
# contains 21,288
#path_2_cat_dr16 = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')
#path_2_cat_dr16 = os.path.join(catalog_dir, '2RXS', 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER_SPIDERSCODEX.fits')

EXI_ML_min_str = sys.argv[1]

if EXI_ML_min_str=='6p5':
	EXI_ML_min=6.5
if EXI_ML_min_str=='10':
	EXI_ML_min=10

prefix='2RXS_ExiML'+EXI_ML_min_str+'_'
#def get_arrays(path_2_cat, EXI_ML_min = EXI_ML_min):
	#hd_rass_i = fits.open(path_2_cat)[1].data
	#targets = (hd_rass_i['NWAY_match_flag']!=2) & (hd_rass_i['RXS_IN_BOSS']==1) & (hd_rass_i['FLAG_SDSSv5b_best']==1) & (hd_rass_i['NWAY_p_any']>=0.01) & (hd_rass_i['NUM_SDSS']>=1)
	##  NWAY_match_flag!=2 &&RXS_IN_BOSS==1 &&FLAG_SDSSv5b_best==1 &&
	#hd_rass = fits.open(path_2_cat)[1].data[targets]
	##print(np.unique(hd_rass['CLASS_BEST']))
	#data = hd_rass
	#ra = data['RXS_RAJ2000']
	#dec = data['RXS_DEJ2000']
	#high_conf = (data['RXS_ExiML'] >= EXI_ML_min )
	#targeted  = (high_conf) & (data['SDSS_FIBER2MAG_i']>=MIN_SDSS_FIBER2MAG_i) & (data['SDSS_FIBER2MAG_i']<=MAX_SDSS_FIBER2MAG_i) & (data['SDSS_MODELMAG_i']>=MIN_SDSS_MODELMAG_i)    
	#observed  = (targeted) & (data['DR16_MEMBER']==1)
	#c1 = (observed) & ((data['Z_BEST']>-0.5) | ((data['DR16_Z']>-0.5) & (data['DR16_Z_ERR']>0)))
	#c2 = (c1) & (data['CONF_BEST']==3)
	#c3 = (c1) & (data['CONF_BEST']==2) & ((data['CLASS_BEST']=='BLAZAR')|(data['CLASS_BEST']=='BLLAC'))
	#c4 = (c1) & (data['DR16_SN_MEDIAN_ALL']>=2) & (data['DR16_ZWARNING']==0 )
	#c5 = (c1) & (data['CONF_BEST']==2) & (data['DR16_ZWARNING']==0 ) & (data['VI_REINSPECT_FLAG'] == 0) & (data['VI_NINSPECTORS']>2)
	#c6 = (c1) & (data['CONF_BEST']==2) & (data['DR16_ZWARNING']==0 ) & (data['VI_AM_RECONCILED']==1)
	#idZ =  (c2) | (c3) | (c4) | (c5) | (c6) 
	#blazars_noZ = (idZ) & ((data['CLASS_BEST']=='BLAZAR')|(data['CLASS_BEST']=='BLLAC')) & (data['CONF_BEST']<3)
	#goodZ = (idZ) & (blazars_noZ == False)
	#redmapper_cluster = (goodZ) & (data['REDMAPPER_Separation']<60) & (abs(data['DR16_Z'] - data['REDMAPPER_Z_LAMBDA'])<0.01)
	#spiders_cluster = (goodZ) & (abs(data['SPIDERSCODEX_SCREEN_CLUZSPEC']-data['Z_BEST'])<0.01)
	#clusters = (redmapper_cluster) | (spiders_cluster)
	#stars = (goodZ) & (data['CLASS_BEST']=='STAR')
	#agnZ = (goodZ) & (clusters==False) &  (stars==False)
	#print(len(ra), len(ra[high_conf]) )
	#return data[(high_conf)], ra[(high_conf)], dec[(high_conf)], targeted[(high_conf)], observed[(high_conf)], goodZ[(high_conf)], idZ[(high_conf)], agnZ[(high_conf)], clusters[(high_conf)], blazars_noZ[(high_conf)], stars[(high_conf)], hd_rass_i

#data, ra, dec, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, all_data = get_arrays(path_2_cat_dr16)
print(data['DR16_Z'], data['DR16Q_Z'], data['Z_BEST'])


qty = 'RXS_SRC_FLUX'
xlabel = r'$\log_{10}(F_X [erg.\, cm^{-2}. s^{-1}])$'
XMIN = -13.
XMAX = -10
xlim_min = -13.  
xlim_max = -10. 
DX=0.5 
logDATA = True

good = (data[qty]>0)
YY = np.log10(data[qty])
bins = np.arange(XMIN, XMAX+DX, DX)
x_bins = (bins[1:]+bins[:-1])/2.

fig_out = os.path.join(figure_dir, prefix+'_xray_hist_dr16.png')
p.figure(1, (5.5,4.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()
p.hist(YY[good ]           , bins=bins, histtype='step',lw=5, label='A')
p.hist(YY[good & targeted] , bins=bins, histtype='step',lw=4, label='T')
p.hist(YY[good & observed] , bins=bins, histtype='step',lw=3, label='O')
p.hist(YY[good & idZ]      , bins=bins, histtype='step',lw=2, label='I')
p.hist(YY[good & goodZ]    , bins=bins, histtype='step',lw=1, label='Z')  
p.xlabel(xlabel)
p.ylabel('N')
p.grid()
p.xlim((xlim_min, xlim_max ))
#p.ylim((0.9,5000))
p.yscale('log')
p.legend(frameon=True, loc=0)
print(fig_out)
p.savefig(fig_out)
p.clf()



qty = 'SDSS_FIBER2MAG_i'
xlabel = 'SDSS FIBER2MAG i' 
XMIN = 10.
XMAX = 30
xlim_min = 10.  
xlim_max = 25. 
DX=0.5 

good = (data[qty]>0)
YY = data[qty]
bins = np.arange(XMIN, XMAX+DX, DX)
x_bins = (bins[1:]+bins[:-1])/2.

fig_out = os.path.join(figure_dir, prefix+'_fiber2magi_hist_dr16.png')
p.figure(1, (5.5,4.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()
p.hist(YY[good ]           , bins=bins, histtype='step',lw=5, label='A')
p.hist(YY[good & targeted] , bins=bins, histtype='step',lw=4, label='T')
p.hist(YY[good & observed] , bins=bins, histtype='step',lw=3, label='O')
p.hist(YY[good & idZ]      , bins=bins, histtype='step',lw=2, label='I')
p.hist(YY[good & goodZ]    , bins=bins, histtype='step',lw=1, label='Z')  
p.xlabel(xlabel)
p.ylabel('N')
p.grid()
p.xlim((xlim_min, xlim_max ))
#p.ylim((0.9,5000))
p.yscale('log')
p.legend(frameon=True, loc=2)
print(fig_out)
p.savefig(fig_out)
p.clf()


qty = 'SDSS_MODELMAG_i'
xlabel = 'SDSS MODELMAG i' 
XMIN = 10.
XMAX = 30
xlim_min = 10.  
xlim_max = 25. 
DX=0.5 

good = (data[qty]>0)
YY = data[qty]
bins = np.arange(XMIN, XMAX+DX, DX)
x_bins = (bins[1:]+bins[:-1])/2.

fig_out = os.path.join(figure_dir, prefix+'_modelmagi_hist_dr16.png')
p.figure(1, (5.5,4.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()
p.hist(YY[good ]           , bins=bins, histtype='step',lw=5, label='A')
p.hist(YY[good & targeted] , bins=bins, histtype='step',lw=4, label='T')
p.hist(YY[good & observed] , bins=bins, histtype='step',lw=3, label='O')
p.hist(YY[good & idZ]      , bins=bins, histtype='step',lw=2, label='I')
p.hist(YY[good & goodZ]    , bins=bins, histtype='step',lw=1, label='Z')  
p.xlabel(xlabel)
p.ylabel('N')
p.grid()
p.xlim((xlim_min, xlim_max ))
#p.ylim((0.9,5000))
p.yscale('log')
p.legend(frameon=True, loc=0)
print(fig_out)
p.savefig(fig_out)
p.clf()

def plot_hist(data=data,  qty='SDSS_FIBER2MAG_i', xlabel='SDSS_FIBER2MAG_i',XMIN = 16,XMAX = 30, xlim_min = 16.8, xlim_max = 22.7, DX= 0.5, logX = False, logDATA = False ):
	good = (data[qty]>0)
	if logDATA:
		YY = np.log10(data[qty])
	else:
		YY = data[qty]
	bins = np.arange(XMIN, XMAX+DX, DX)# (np.max(YY)-np.min(YY))/10.)
	x_bins = (bins[1:]+bins[:-1])/2.
	#
	out_all = np.histogram( YY[good & observed] , bins=bins)[0]
	out_IDZ = np.histogram( YY[good & idZ]      , bins=bins)[0]
	out_BLA = np.histogram( YY[good & blazars_noZ], bins=bins)[0]
	out_GOO = np.histogram( YY[good & goodZ]    , bins=bins)[0]
	out_AGN = np.histogram( YY[good & agnZ]     , bins=bins)[0]
	out_CLU = np.histogram( YY[good & clusters] , bins=bins)[0]
	out_STA = np.histogram( YY[good & stars]    , bins=bins)[0]
	#
	fig_out = os.path.join(figure_dir, prefix+'_'+qty+'_hist_dr16.png')
	p.figure(1, (5.5,4.5))
	p.axes([0.15, 0.15, 0.8, 0.77])
	p.tight_layout()
	p.hist(YY[good & observed] , bins=bins, histtype='step',lw=4, label='O')
	p.hist(YY[good & idZ]      , bins=bins, histtype='step',lw=3, label='I')
	#p.hist(YY[good & blazars_noZ], bins=bins, histtype='step',lw=3, label='blzr')
	p.hist(YY[good & goodZ]    , bins=bins, histtype='step',lw=2, label='Z')  
	#p.hist(YY[good & agnZ]     , bins=bins, histtype='step',lw=3, label='AGN')
	#p.hist(YY[good & clusters] , bins=bins, histtype='step',lw=3, label='cluster')
	#p.hist(YY[good & stars]    , bins=bins, histtype='step',lw=3, label='star') 
	p.xlabel(xlabel)
	p.ylabel('N')
	p.grid()
	p.xlim((xlim_min, xlim_max ))
	#p.ylim((0.9,5000))
	#p.yscale('log')
	p.legend(frameon=True, loc=0)
	print(fig_out)
	p.savefig(fig_out)
	p.clf()

	# RATIO with observed
	fig_out = os.path.join(figure_dir, prefix+'_'+qty+'_ratio_dr16.png')
	p.figure(1, (5.5,4.5))
	p.axes([0.15, 0.15, 0.8, 0.77])
	p.tight_layout()
	# blazar noZ
	mid_line = out_BLA*1./out_all
	err_line = out_BLA**(-0.5)
	p.errorbar(x_bins, y = mid_line, xerr=DX/2., label='Blzr/O', fmt='none', lw=3)
	# identified
	mid_line = out_IDZ*1./out_all
	err_line = out_IDZ**(-0.5)
	p.errorbar(x_bins, y = mid_line, xerr=DX/2., label='I/O', fmt='none', lw=3)
	# redshift
	mid_line = out_GOO*1./out_all
	err_line = out_GOO**(-0.5)
	p.errorbar(x_bins, y = mid_line, xerr=DX/2., label='Z/O', fmt='none', lw=3)
	#
	p.xlabel(xlabel)
	p.ylabel('fraction')
	p.grid()
	p.xlim((xlim_min, xlim_max ))
	p.ylim((-0.01,1.01))
	#p.yscale('log')
	p.legend(frameon=True, loc=6)
	print(fig_out)
	p.savefig(fig_out)
	p.clf()

	# RATIO with goodZ
	fig_out = os.path.join(figure_dir, prefix+'_'+qty+'_ratioZ_dr16.png')
	p.figure(1, (5.5,4.5))
	p.axes([0.15, 0.15, 0.8, 0.77])
	p.tight_layout()
	# blazar noZ
	mid_line = out_AGN*1./out_GOO
	err_line = out_AGN**(-0.5)
	p.errorbar(x_bins, y = mid_line, xerr=DX/2., label='AGN/Z', fmt='none', lw=3)
	# identified
	mid_line = out_CLU*1./out_GOO
	err_line = out_CLU**(-0.5)
	p.errorbar(x_bins, y = mid_line, xerr=DX/2., label='cluster/Z', fmt='none', lw=3)
	# redshift
	mid_line = out_STA*1./out_GOO
	err_line = out_STA**(-0.5)
	p.errorbar(x_bins, y = mid_line, xerr=DX/2., label='star/Z', fmt='none', lw=3)
	#
	p.xlabel(xlabel)
	p.ylabel('fraction')
	p.grid()
	p.xlim((xlim_min, xlim_max ))
	p.ylim((-0.01,1.01))
	if logX : #= False
		p.xscale('log')
	#p.yscale('log')
	p.legend(frameon=True, loc=7)
	print(fig_out)
	p.savefig(fig_out)
	p.clf()


plot_hist(data=data, 
			  qty='SDSS_FIBER2MAG_i',
			  xlabel = 'SDSS FIBER2MAG i',
			  XMIN = 16,
			  XMAX = 30, 
			  xlim_min = 16.8, 
			  xlim_max = 22.7, 
			  DX= 0.5 )

plot_hist(data=data, 
			  qty='SDSS_MODELMAG_i',
			  xlabel = 'SDSS MODELMAG i',
			  XMIN = 16,
			  XMAX = 30, 
			  xlim_min = 15.8, 
			  xlim_max = 22.2, 
			  DX= 0.5 )

plot_hist(data=data, 
			  qty='SDSS_MODELMAG_r', 
			  xlabel = 'SDSS MODELMAG r',
			  XMIN = 16,
			  XMAX = 30, 
			  xlim_min = 15.8, 
			  xlim_max = 22.2, 
			  DX= 0.5 )

plot_hist(data=data, 
			  qty='DR16_SN_MEDIAN_ALL', 
			  xlabel = 'MEDIAN(S/N)',
			  XMIN = -1,
			  XMAX = 50, 
			  xlim_min = 0.4,  
			  xlim_max = 20, 
			  DX= 1., 
			  logX = True)

plot_hist(data=data, 
		qty = 'RXS_SRC_FLUX',
		xlabel = r'$\log_{10}(F_X [erg.\, cm^{-2}. s^{-1}])$',
		XMIN = -13.,
		XMAX = -10,
		xlim_min = -13.,  
		xlim_max = -10., 
		DX=0.5, logDATA = True)

