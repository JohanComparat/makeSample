import os, sys
import numpy as np
import matplotlib.pyplot as p

from scipy.stats import kstest

out_dir = '/data36s/comparat/AGN_clustering/angular_clustering/'

def compute_clustering(name, out_dir):
	out_name = os.path.join(out_dir , name )
	f=open(out_name+'.ini','w')
	f.write('data_filename= '+out_name+'.data \n')    
	f.write('random_filename= '+out_name+'.random \n')    
	f.write('input_format= 2 \n')
	f.write('output_filename= '+out_name+'.wtheta \n')
	f.write('corr_type= angular \n')
	f.write('omega_M= 0.307 \n')
	f.write('omega_L= 0.693 \n')
	f.write('w= -1 \n')
	f.write('log_bin= 10 \n')
	f.write('dim1_min_logbin= 0.001 \n')
	f.write('dim1_max= 1. \n')
	f.write('dim1_nbin= 60 \n')
	f.write('dim2_max= 160. \n')
	f.write('dim2_nbin= 40 \n')
	f.write('dim3_min= 0.00 \n')
	f.write('dim3_max= 3. \n')
	f.write('dim3_nbin= 1 \n')
	f.write('use_pm= 0 \n')
	f.write('n_pix_sph= 1000000 \n')
	f.close()
	os.system("~/darksim/software/CUTE/CUTE/CUTE "+out_name+'.ini')



def plot_results(name, out_dir):
	out_name = os.path.join(out_dir , name )

	fig = p.figure(3, (9,12))
	# wedge plot
	fig.add_subplot(3,2,1)

	ra_D, dec_D, z_D, w = np.loadtxt(out_name + '.data', unpack=True)
	p.plot(ra_D, dec_D, marker=',', color='b', alpha=0.5, rasterized = True, ls='None', label=str(len(ra_D)))
	p.xlabel('R.A. [deg]')
	p.ylabel('Dec. [deg]')
	p.legend(frameon=False, loc=0)
	p.title('data')

	fig.add_subplot(3,2,2)
	ra_R, dec_R, z, w = np.loadtxt(out_name + '.random', unpack=True)
	p.plot(ra_R, dec_R, marker=',', color='b', alpha=0.025, rasterized = True, ls='None', label=str(len(ra_R)))
	p.xlabel('R.A. [deg]')
	p.ylabel('Dec. [deg]')
	p.legend(frameon=False, loc=0)
	p.title('randoms')

	fig.add_subplot(3,2,3)
	# KS-test RA
	nn_R,bb,pp=p.hist(ra_R, bins=100, histtype='step', cumulative=True, normed=True, lw=4)
	nn_D,bb,pp=p.hist(ra_D, bins=bb, histtype='step', cumulative=True, normed=True, lw=2)
	KS_out = np.max(nn_R-nn_D)
	p.title(str(KS_out))
	p.xlabel('R.A. [deg]')
	p.ylabel('cdf')
	p.grid()

	fig.add_subplot(3,2,4)
	# KS-test DEC
	nn_R,bb,pp=p.hist(dec_R, bins=100, histtype='step', cumulative=True, normed=True, lw=4)
	nn_D,bb,pp=p.hist(dec_D, bins=bb, histtype='step', cumulative=True, normed=True, lw=2)
	KS_out = np.max(nn_R-nn_D)
	p.title(str(KS_out))
	p.ylabel('Dec. [deg]')
	p.ylabel('cdf')
	p.grid()

	fig.add_subplot(3,2,5)
	# nz
	DATA = np.loadtxt(out_name + '.wtheta', unpack=True)
	p.errorbar(DATA[0], DATA[1], yerr=DATA[1]*DATA[2]**(-0.5), rasterized = True)
	p.plot(DATA[0], 0.1*DATA[0]**-0.8)
	p.xlabel('theta [deg]')
	p.ylabel('theta x w(theta) ')
	p.xscale('log')
	p.yscale('log')
	p.grid()

	fig.add_subplot(3,2,6)
	# nz
	p.hist(z_D[z_D>0], bins=np.arange(0,2.6,0.1), rasterized = True)
	p.xlabel('redshift')
	p.ylabel('Counts / dz=0.1')
	#p.xscale('log')
	p.yscale('log')
	p.title(str(np.median(z_D[z_D>0])))
	p.grid()
	
	p.tight_layout()
	p.savefig(out_name+".png")
	p.clf()

#compute_clustering('XMMSL2_AllWISE_catalog_paper_2017JUN09_X_GAL', out_dir)
#plot_results('XMMSL2_AllWISE_catalog_paper_2017JUN09_X_GAL', out_dir)

#compute_clustering('2RXS_AllWISE_catalog_paper_2017May26_X_GAL', out_dir)
#plot_results('2RXS_AllWISE_catalog_paper_2017May26_X_GAL', out_dir)

#compute_clustering('2RXS_AllWISE_catalog_paper_2017May26_X_GAL_noStars', out_dir)
#plot_results('2RXS_AllWISE_catalog_paper_2017May26_X_GAL_noStars', out_dir)

#compute_clustering('2RXS_AllWISE_catalog_paper_2017May26_X_GAL_noStars_rtlim_gt_0_005', out_dir)
#plot_results('2RXS_AllWISE_catalog_paper_2017May26_X_GAL_noStars_rtlim_gt_0_005', out_dir)

#compute_clustering('2RXS_AllWISE_catalog_paper_2017May26_X_GAL_noStars_rtlim_gt_0_015', out_dir)
#plot_results('2RXS_AllWISE_catalog_paper_2017May26_X_GAL_noStars_rtlim_gt_0_015', out_dir)

#compute_clustering('2RXSALLWISE_XGAL_ratelim004_gaia12', out_dir)
plot_results('2RXSALLWISE_XGAL_ratelim004_gaia12', out_dir)
os.system('cp /data36s/comparat/AGN_clustering/angular_clustering/*.png /home/comparat/wwwDir/stuff/')

#ls /data36s/comparat/AGN_clustering/angular_clustering/*.random

all_files = np.array([
'2RXSALLWISE_XGAL_ratelim_005_gaia12'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_10'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_11'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_12'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_13'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_14'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_15'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_16'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_17'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_7'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_8'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_mask_gaia_g_lt_9'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_noStars'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_noStars_rtlim_gt_0_005'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL_noStars_rtlim_gt_0_015'
#,'2RXS_AllWISE_catalog_paper_2017May26_X_GAL'
#,'XMMSL2_AllWISE_catalog_paper_2017JUN09_X_GAL'
])


for ffi in all_files:
	#compute_clustering(ffi, out_dir)
	plot_results(ffi, out_dir)

os.system('cp /data36s/comparat/AGN_clustering/angular_clustering/*.png /home/comparat/wwwDir/stuff/')
