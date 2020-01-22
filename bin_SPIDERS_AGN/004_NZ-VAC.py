"""
python plot_NZ-VAC.py MD04
python plot_NZ-VAC.py MD10

"""
from catalog_lib import *

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

env = "MD10" # sys.argv[1]

area_eboss_dr16 = 5128.

#sig_fx = 0.3 #1e-13
#sig_fx = 0.2 #1e-13
#sig_fx = 0.1 #sig_fx = 2e-13
sig_fx = float(sys.argv[1] )

MAG_MIN_LIM = 15
#MAG_LIM = 22.5
MAG_LIM = float(sys.argv[2] )

#EXI_ML_min_str = sys.argv[2]
#if EXI_ML_min_str=='6p5':
#EXI_ML_min=6.5
#if EXI_ML_min_str=='10':
#EXI_ML_min=10
EXI_ML_min = float(sys.argv[3] )

redshift = data['Z_BEST']
redshift[redshift<0]=0
#Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2]

dz=0.05
z_min=0.0005 
z_max = np.max(redshift[agnZ]) # np.max([np.max(redshift),np.max(Z_XMMSL2)])
zs=np.arange(z_min, z_max + dz, dz)
xzs=(zs[1:]+zs[:-1])/2.

def get_z_arrs_GAL(path_2_mock):
	print(path_2_mock)
	hd_mock_a = fits.open(path_2_mock)[1].data
	lx = hd_mock_a['galaxy_LX_hard'] 
	dl_cm = dL_interpolation(hd_mock_a['redshift_R'])
	fx_true_a = 10*10**lx / (4. * np.pi * dl_cm**2)
	filtermock = (hd_mock_a['redshift_R']<z_max)&(fx_true_a>2e-15)
	hd_mock = hd_mock_a[filtermock]
	fx_true = fx_true_a[filtermock]
	fx_logerr = norm.rvs(loc = 0, scale=sig_fx, size=len(fx_true)  )
	fx = 10**( np.log10(fx_true) + fx_logerr )
	Z_mock_1em125_a = np.histogram( hd_mock['redshift_R'][(fx>10**(-12.5)) ] , bins=zs)[0]
	Z_mock_1em125_1 = np.histogram( hd_mock['redshift_R'][(fx>10**(-12.6)) ] , bins=zs)[0]
	Z_mock_1em125_2 = np.histogram( hd_mock['redshift_R'][(fx>10**(-12.7)) ] , bins=zs)[0]
	Z_mock_1em13_a  = np.histogram( hd_mock['redshift_R'][(fx>10**(-12.8)) ] , bins=zs)[0]
	Z_mock_1em13_1  = np.histogram( hd_mock['redshift_R'][(fx>10**(-13.0)) ] , bins=zs)[0]
	Z_mock_1em13_2  = np.histogram( hd_mock['redshift_R'][(fx>10**(-13.1)) ] , bins=zs)[0]
	return Z_mock_1em125_a, Z_mock_1em125_1, Z_mock_1em125_2, Z_mock_1em13_a, Z_mock_1em13_1, Z_mock_1em13_2


sdss_depth_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'sdss_depth')
depth_sdss_r = os.path.join(
    sdss_depth_dir,
    'sdss_dr8_nodered_nside2048_r_model_10sigma.fits.gz')

NSIDE = 2048
sdss_r_depth = hp.fitsfunc.read_map(depth_sdss_r)
median_depth = np.median(sdss_r_depth[sdss_r_depth>0])

# conver to flux
def get_sigma_m(depth_10_sigma, sdss_temp):
	limiting_flux = 2.5 * 10**(-0.4 * (depth_10_sigma)) / np.log(10)/10.
	observed_flux = 2.5 * 10**(-0.4 * (sdss_temp)) / np.log(10)
	return 2.5 * limiting_flux / observed_flux / np.log(10)


def get_z_arrs(path_2_mock, output_data=False):
	print(path_2_mock)
	hd_mock_a = fits.open(path_2_mock)[1].data
	fx_true_a = hd_mock_a['agn_FX_soft'] 
	hd_mock = hd_mock_a[fx_true_a>2e-15]
	fx_true = hd_mock['agn_FX_soft'] 
	redshift = hd_mock['redshift_R']
	LX_soft = hd_mock['AGN_LX_soft']
	N_obj = len(fx_true)
	fx_logerr = norm.rvs(loc = 0, scale=sig_fx, size=N_obj  )
	fx = 10**( np.log10(fx_true) + fx_logerr )
	sdss_r_temp = hd_mock['AGN_SDSS_r_magnitude']
	# retrieves the depths
	#pixels = healpy.ang2pix(NSIDE,hd_mock['ra'],hd_mock['dec'], nest=False, lonlat=True)
	depth_10_sigma_r = np.ones_like(sdss_r_temp) * median_depth #sdss_r_depth[pixels]
	# computes the magnitudes
	# their uncertainty
	sdss_r_err = get_sigma_m(depth_10_sigma_r, sdss_r_temp)
	sdss_r_err[sdss_r_err == np.inf] = 0
	rds = norm.rvs(loc=0, scale=1, size=N_obj)
	mag = sdss_r_temp + rds * sdss_r_err
	AGN_type = hd_mock['AGN_type']
	t1 = (AGN_type==11)|(AGN_type==12) 
	t2 = (AGN_type==21)|(AGN_type==22) 
	Z_mock_1em125_a  = np.histogram( redshift[(fx>10**(-12.6+redshift/5.)) & (mag>15) & (mag<MAG_LIM)] , bins=zs)[0]
	Z_mock_1em125_1  = np.histogram( redshift[(fx>10**(-12.6+redshift/5.)) & (mag>15) & (mag<MAG_LIM) & (t1)] , bins=zs)[0]
	Z_mock_1em127_a  = np.histogram( redshift[(fx>10**(-12.7+redshift/5.)) & (mag>15) & (mag<MAG_LIM)] , bins=zs)[0]
	Z_mock_1em127_1  = np.histogram( redshift[(fx>10**(-12.7+redshift/5.)) & (mag>15) & (mag<MAG_LIM) & (t1)] , bins=zs)[0]
	Z_mock_1em129_a  = np.histogram( redshift[(fx>10**(-12.8+redshift/5.)) & (mag>15) & (mag<MAG_LIM)] , bins=zs)[0]
	Z_mock_1em129_1  = np.histogram( redshift[(fx>10**(-12.8+redshift/5.)) & (mag>15) & (mag<MAG_LIM) & (t1)] , bins=zs)[0]
	if output_data:
		return Z_mock_1em125_a, Z_mock_1em125_1, Z_mock_1em127_a, Z_mock_1em127_1, Z_mock_1em129_a, Z_mock_1em129_1, fx, mag, redshift, LX_soft
	else:
		return Z_mock_1em125_a, Z_mock_1em125_1, Z_mock_1em127_a, Z_mock_1em127_1, Z_mock_1em129_a, Z_mock_1em129_1


# AGN NZ from the mock catalogue
path_2_mock = os.path.join(os.environ[env],'SPIDERS_AGN_all.fit')
#sim_list5 = np.array(glob.glob(os.path.join(os.environ[env],'cat_AGN_all/00035?.fit')))
#sim_list6=np.array(glob.glob(os.path.join(os.environ[env],'cat_AGN_all/00045?.fit')))
#sim_list = np.hstack((sim_list5,sim_list6))
#simulation_area = healpy.nside2pixarea(8, degrees=True) * len(sim_list)
#ZS_data = np.sum(np.array([get_z_arrs(p2s) for p2s in sim_list]), axis=0)/simulation_area
#path_2_mock = sim_list[0]
#Z_mock_1em125_a, Z_mock_1em125_1, Z_mock_1em127_a, Z_mock_1em127_1, Z_mock_1em129_a, Z_mock_1em129_1, mock_fx, mock_mag, mock_redshift, mock_LX_soft = get_z_arrs(path_2_mock, output_data=True)

print(path_2_mock)
hd_mock_a = fits.open(path_2_mock)[1].data
fx_true_a = hd_mock_a['agn_FX_soft'] 
simulation_area = 30575. #129600/np.pi
hd_mock = hd_mock_a[ ( fx_true_a > 2e-15 ) & ( abs(hd_mock_a['g_lat']) > 15 ) ]
fx_true = hd_mock['agn_FX_soft'] 
redshift_R = hd_mock['redshift_R']
LX_soft = hd_mock['AGN_LX_soft']
N_obj = len(fx_true)
fx_logerr = norm.rvs(loc = 0, scale=sig_fx, size=N_obj  )
fx = 10**( np.log10(fx_true) + fx_logerr )
sdss_r_temp = hd_mock['AGN_SDSS_r_magnitude']
# retrieves the depths
#pixels = healpy.ang2pix(NSIDE,hd_mock['ra'],hd_mock['dec'], nest=False, lonlat=True)
depth_10_sigma_r = np.ones_like(sdss_r_temp) * median_depth #sdss_r_depth[pixels]
# computes the magnitudes
# their uncertainty
sdss_r_err = get_sigma_m(depth_10_sigma_r, sdss_r_temp)
sdss_r_err[sdss_r_err == np.inf] = 0
rds = norm.rvs(loc=0, scale=1, size=N_obj)
mag = sdss_r_temp + rds * sdss_r_err
AGN_type = hd_mock['AGN_type']
t1 = (AGN_type==11)|(AGN_type==12) 
t2 = (AGN_type==21)|(AGN_type==22) 
Z_mock_low_a  = np.histogram( redshift_R[(fx>10**(-12.95+redshift_R/10.)) & (mag>MAG_MIN_LIM) & (mag<MAG_LIM)] , bins=zs)[0]/simulation_area
Z_mock_low_1  = np.histogram( redshift_R[(fx>10**(-12.95+redshift_R/10.)) & (mag>MAG_MIN_LIM) & (mag<MAG_LIM) & (t1)] , bins=zs)[0]/simulation_area
Z_mock_low_2  = np.histogram( redshift_R[(fx>10**(-12.95+redshift_R/10.)) & (mag>MAG_MIN_LIM) & (mag<MAG_LIM) & (t2)] , bins=zs)[0]/simulation_area

Z_mock_mid_a  = np.histogram( redshift_R[(fx>10**(-13.0+redshift_R/10.)) & (mag>MAG_MIN_LIM) & (mag<MAG_LIM)] , bins=zs)[0]/simulation_area
Z_mock_mid_1  = np.histogram( redshift_R[(fx>10**(-13.0+redshift_R/10.)) & (mag>MAG_MIN_LIM) & (mag<MAG_LIM) & (t1)] , bins=zs)[0]/simulation_area
Z_mock_mid_2  = np.histogram( redshift_R[(fx>10**(-13.0+redshift_R/10.)) & (mag>MAG_MIN_LIM) & (mag<MAG_LIM) & (t2)] , bins=zs)[0]/simulation_area

Z_mock_up_a  = np.histogram( redshift_R[(fx>10**(-13.05+redshift_R/10.)) & (mag>MAG_MIN_LIM) & (mag<MAG_LIM)] , bins=zs)[0]/simulation_area
Z_mock_up_1  = np.histogram( redshift_R[(fx>10**(-13.05+redshift_R/10.)) & (mag>MAG_MIN_LIM) & (mag<MAG_LIM) & (t1)] , bins=zs)[0]/simulation_area
Z_mock_up_2  = np.histogram( redshift_R[(fx>10**(-13.05+redshift_R/10.)) & (mag>MAG_MIN_LIM) & (mag<MAG_LIM) & (t2)] , bins=zs)[0]/simulation_area

#if output_data:
#return Z_mock_1em125_a, Z_mock_1em125_1, Z_mock_1em127_a, Z_mock_1em127_1, Z_mock_1em129_a, Z_mock_1em129_1, fx, mag, redshift_R, LX_soft
#else:
#ZS_data = Z_mock_1em125_a, Z_mock_1em125_1, Z_mock_1em127_a, Z_mock_1em127_1, Z_mock_1em129_a, Z_mock_1em129_1


# GALAXY NZ from the mock catalogue
#sim_list_gal=np.array(glob.glob(os.path.join(os.environ[env],'cat_GALAXY_all/00035?.fit')))[:4]
#simulation_area = healpy.nside2pixarea(8, degrees=True) * len(sim_list_gal)
#ZS_data_GAL = np.sum(np.array([get_z_arrs_GAL(p2s) for p2s in sim_list_gal]), axis=0)/simulation_area
#path_2_mock = sim_list_gal[0]

# REDSHIFT HISTOGRAM
p.figure(0, (10,5))
p.axes([0.12, 0.12, 0.78, 0.78])
# AGN mocks
p.fill_between(xzs, y1=Z_mock_low_a, y2=Z_mock_up_a, color='b', alpha=0.1, label='C19 mock')
#p.fill_between(xzs, y1=Z_mock_low_1, y2=Z_mock_up_1, color='orange', alpha=0.1, label='t1')
#p.fill_between(xzs, y1=Z_mock_low_2, y2=Z_mock_up_2, color='orange', alpha=0.1, label='mock t2')
#p.plot(xzs, ZS_data[0], ls='solid', lw=1,  label='A, FX>-12.95')
#p.plot(xzs, ZS_data[1], ls='dashed', lw=1, label='A, FX>-12.95 t1')
#p.plot(xzs, ZS_data[2], ls='solid', lw=1,  label='A, FX>-13.0')
#p.plot(xzs, ZS_data[3], ls='dashed', lw=1, label='A, FX>-13.0 t1')
#p.plot(xzs, ZS_data[4], ls='solid', lw=1,  label='A, FX>-13.05')
#p.plot(xzs, ZS_data[5], ls='dashed', lw=1, label='A, FX>-13.05 t1')
# GALAXY mock
#p.plot(xzs, ZS_data_GAL[0], ls='solid', lw=2,  label='G, FX>-12.5')
#p.plot(xzs, ZS_data_GAL[1], ls='dashed', lw=2, label='G, FX>-12.6')
#p.plot(xzs, ZS_data_GAL[2], ls='dotted', lw=2, label='G, FX>-12.7')
#p.plot(xzs, ZS_data_GAL[3], ls='solid', lw=2,  label='G, FX>-12.8')
#p.plot(xzs, ZS_data_GAL[4], ls='dashed', lw=2, label='G, FX>-13.0')
#p.plot(xzs, ZS_data_GAL[5], ls='dotted', lw=2, label='G, FX>-13.1')
# 2RXS SPIDERS DATA
# All AGN
NN_2rxs = np.histogram(redshift[agnZ], bins=zs)[0]
density = NN_2rxs/area_eboss_dr16#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs**(-0.5), xerr=dz/2., label='AGN 2RXS', fmt='s', mfc='none')

NN_clusters = np.histogram(redshift[clusters], bins=zs)[0]
density = NN_clusters/area_eboss_dr16#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_clusters**(-0.5), xerr=dz/2., label='Clusters 2RXS', fmt='s', mfc='none')

# REDSHIFT HISTOGRAM saved as ascii
#np.savetxt(os.path.join( figure_dir, env+"_histogram_redshift_"+str(sig_fx)+'_'+str(MAG_LIM)+'_'+str(EXI_ML_min)+".txt"), np.transpose([zs[:-1], zs[1:], NN_2rxs]), delimiter=" & ", newline=" \\\\ \n", fmt='%10.1f')

dz = 0.1
zs = np.arange(z_min, z_max + dz, dz)
xzs = (zs[1:]+zs[:-1])/2.

NN_2rxs = np.histogram(redshift[agnZ], bins=zs)[0]
density = NN_2rxs/area_eboss_dr16#/np.max(NN)

NN_clusters = np.histogram(redshift[clusters], bins=zs)[0]
density = NN_clusters/area_eboss_dr16#/np.max(NN)

np.savetxt(os.path.join( figure_dir, env+"_histogram_redshift_"+str(sig_fx)+'_'+str(MAG_LIM)+'_'+str(EXI_ML_min)+".txt"), np.transpose([zs[:-1], zs[1:], NN_2rxs]), delimiter=" & ", newline=" \\\\ \n", fmt='%10.1f')

print('NN_2rxs SUM', np.sum(NN_2rxs))

dz=0.05
z_min=0.0005 
z_max = np.max(redshift[agnZ]) # np.max([np.max(redshift),np.max(Z_XMMSL2)])
zs=np.arange(z_min, z_max + dz, dz)
xzs=(zs[1:]+zs[:-1])/2.

data, high_conf, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, bl, nl, class_best = get_arrays_xmmsl2(path_2_cat_xmmsl = path_2_cat_xmmsl, EXI_ML_min = EXI_ML_min)


redshift = data['Z_BEST']

NN_2rxs = np.histogram(redshift[agnZ], bins=zs)[0]
density = NN_2rxs/area_eboss_dr16#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs**(-0.5), xerr=dz/2., label='AGN XMMSL2', fmt='o', mfc='none')

NN_clusters = np.histogram(redshift[clusters], bins=zs)[0]
density = NN_clusters/area_eboss_dr16#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_clusters**(-0.5), xerr=dz/2., label='Clusters XMMSL2', fmt='o', mfc='none')

print('NN_2rxs SUM', np.sum(NN_2rxs))


## Type 1
#NN_2rxs_t1 = np.histogram(redshift[bl], bins=zs)[0]
#density = NN_2rxs_t1/area_eboss_dr16#/np.max(NN)
#p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs_t1**(-0.5), xerr=dz/2., label='2RXS t1', fmt='none')
# Type 2
#NN_2rxs_t2 = np.histogram(redshift[nl], bins=zs)[0]
#density = NN_2rxs_t2/area_eboss_dr16#/np.max(NN)
#p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs_t2**(-0.5), xerr=dz/2., label='2RXS t2', fmt='none')
# labels
p.xlabel('redshift')
p.ylabel('N/deg2, dz='+str(dz))
p.yscale('log')
p.xlim((0,3.2))
p.grid()
p.legend(frameon=False, loc=0)
#p.title(str(sig_fx)+'_'+str(MAG_LIM)+'_'+str(EXI_ML_min) )
fig_out = os.path.join( figure_dir, env+"_histogram_redshift_"+str(sig_fx)+'_'+str(MAG_LIM)+'_'+str(EXI_ML_min)+".png")
print(fig_out)
p.savefig(fig_out)
p.clf()




#p.figure(0, (10,5))
#p.axes([0.12, 0.12, 0.78, 0.78])
## AGN mocks
 ##mock_fx, mock_mag, mock_redshift
#dls = np.log10(dL_interpolation(redshift[agnZ]))
#p.plot(redshift[agnZ], data['RXS_LOGLX'][agnZ] -2*dls-np.log10(4*np.pi), 'r+', alpha=0.1)
## 2RXS data
#mock_selection = (fx>10**(-13. + redshift_R/10)) & (mag>13) & (mag<MAG_LIM)
#p.plot(redshift_R[mock_selection], np.log10(fx_true[mock_selection]), 'k,')
##
#p.plot(redshift[agnZ], -12.8 + redshift[agnZ]/10, 'g,')
#p.xlabel('redshift')
#p.ylabel('X-ray flux')
##p.yscale('log')
#p.xlim((0,3))
#p.grid()
#p.legend(frameon=False, loc=0)
#p.title(str(sig_fx)+'_'+str(MAG_LIM)+'_'+str(EXI_ML_min) )
#p.savefig(os.path.join( figure_dir, env+"_LX_vs_redshift_"+str(sig_fx)+'_'+str(MAG_LIM)+'_'+str(EXI_ML_min)+".png"))
#p.clf()















sys.exit()

p.figure(0, (10,5))
p.axes([0.12, 0.12, 0.78, 0.78])

p.plot(xzs, ZS_data[3], ls='solid', lw=2, label='m, FX>-13')
p.plot(xzs, ZS_data[4], ls='dashed', lw=2, label='m, FX>-13 t1')
p.plot(xzs, ZS_data[5], ls='dotted', lw=2, label='m, FX>-13 t2')
#
Z_RASS = redshift[goodZ]
#
NN_2rxs = np.histogram(Z_RASS, bins=zs)[0]
density = NN_2rxs/area_eboss_dr16#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs**(-0.5), xerr=dz/2., label='2RXS', fmt='none')

# Type 1
Z_RASS = redshift[goodZ & bl]
#
NN_2rxs_t1 = np.histogram(Z_RASS, bins=zs)[0]
density = NN_2rxs_t1/area_eboss_dr16#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs_t1**(-0.5), xerr=dz/2., label='2RXS t1', fmt='none')

# Type 2
Z_RASS = redshift[goodZ & nl]
#
NN_2rxs_t2 = np.histogram(Z_RASS, bins=zs)[0]
density = NN_2rxs_t2/area_eboss_dr16#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_2rxs_t2**(-0.5), xerr=dz/2., label='2RXS t2', fmt='none')

p.xlabel('redshift')
p.ylabel('N/deg2  dz=0.2 ')
p.yscale('log')
p.xlim((0,4))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( figure_dir, env+"_histogram_redshift_2RXS_"+str(sig_fx)+".png"))
p.clf()

sys.exit()

p.figure(0, (10,5))
p.axes([0.12, 0.12, 0.78, 0.78])

p.plot(xzs, ZS_data[0], ls='solid', lw=2, label='m, FX>-12.5')
p.plot(xzs, ZS_data[1], ls='dashed', lw=2, label='m, FX>-12.5 t1')
p.plot(xzs, ZS_data[2], ls='dotted', lw=2, label='m, FX>-12.5 t2')
#
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2]
#
NN_xmmxxl = np.histogram(Z_XMMSL2, bins=zs)[0]
density = NN_xmmxxl/area_eboss_dr16#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_xmmxxl**(-0.5), xerr=dz/2., label='XMMSL2', fmt='none')

# Type 1
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2 & bl2]
#
NN_xmmxxl_t1 = np.histogram(Z_XMMSL2, bins=zs)[0]
density = NN_xmmxxl_t1/area_eboss_dr16#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_xmmxxl_t1**(-0.5), xerr=dz/2., label='XMMSL2 t1', fmt='none')

# Type 2
Z_XMMSL2 = hd_xmmsl['Z_BEST'][goodZ2 & nl2]
#
NN_xmmxxl_t2 = np.histogram(Z_XMMSL2, bins=zs)[0]
density = NN_xmmxxl_t2/area_eboss_dr16#/np.max(NN)
p.errorbar((zs[1:]+zs[:-1])/2., density, yerr=density*NN_xmmxxl_t2**(-0.5), xerr=dz/2., label='XMMSL2 t2', fmt='none')

p.xlabel('redshift')
p.ylabel('N/deg2  dz=0.2 ')
p.yscale('log')
p.xlim((0,4))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( figure_dir, env+"_histogram_redshift_XMMSL2.png"))
p.clf()

