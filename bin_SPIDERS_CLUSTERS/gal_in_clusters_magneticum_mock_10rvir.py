"""
Omega_M=0.272
Omega_B=0.0456      (= 16.8 %)
Omega_L=0.728
h=0.704
n=0.963
sigma_8=0.809

---------------------- Lightcone geometry -- geometry.dat ---------------------------

This is a 1/8th of a sky-wedge, reaching out to z~1.25. The wedge consist of 23
independent slices between z=0.0 and z=1.25, the geometrical configuration can be
found in wmap.geometry.dat. For z up to 0.21372931 Box2/hr was used, for redshift above
Box2b/hr was used. Note that the z=0.21372931 slice contains contributions from Box2/hr
(snap 116) and Box2b/hr (snap 031). Every Box is duplicated five times along each direction
to fit in the wedges, am display of the geometry can be found in sketch_outer_part.jpg. 

The file wmap.geometry.dat contains:
snapNr:   Snap Number
z:        redshift of the simulation output used for the slice
r0:       proper distance to the start of the wedge
r1:       proper distance to the end of the wedge


Clusters
========

ihal:        HaloID
l,b:         position on the sky of the cluster (in degree)
r:           proper distance in MPc
z_true:      real redshift of the object (use this for distances)
z_obs:       redshift as spectroscopicly observed (use this as observed redshift)
Mvir:        Mvir (spherical tophat model) in Msol/h
Rvir:        Rvir in kpc/h.
M500:        M500 (critical) in Msol/h
R500:        R500 (critical) in kpc/h.
Mstar500:    Stellar mass within R500 (critical) in Msol/h
Mgas500:     Gas mass within R500 (critical) in Msol/h
T500:        Mass weighted temperature within R500 (critical) in keV
Lx500:       Bolometric X-Ray luminosity within R500 (critical) in 1e44 erg/s
Ysz500:      ComptonY within R500 (critical)
M200:        M200 (critical) in Msol/h
R200         R200 (critical) in kpc/h

Galaxies
========

isub l b rr vmax z_true z_obs Mstar sfr u  V  g  r  i  z  Y  J  H  K  L  M Age Z
isub:        SubHaloID
l,b:         position on the sky of the cluster (in degree)
rr:          proper distance in MPc
vmax:        maximum circular velocity of the dark matter halo.
z_true:      real redshift of the object (use this for distances)
z_obs:       redshift as spectroscopicly observed (use this as observed redshift)
Mstar:       stellar mass in Msol/h
sfr:         starformation rate in Msol/year

u,V,g,r,i,z,Y,J,H,K,L,M: 
magnitudes based on BS2007 with a underlying Charbrier IMF in VEGA system. Given
are the observer frame where we do take attenuation into account (we here follow
the description by DeLucia). You also need to use luminosity distance to convert to
observed magnitudes.

Age:         Mean age of the stellar population (expressed as the expansion factor)
Z:           Mean metalicity of the stellar component

"""
#######################
# PACKAGES
#######################

import os 
import numpy as n

import pickle 
import astropy.io.fits as fits

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

from scipy.interpolate import interp1d

import time
t0 = time.time()
print('start', t0)
#######################
# COSMO
#######################

import astropy.cosmology as co
import astropy.units as u
cosmo = co.WMAP7

dC_2_z = interp1d(cosmo.comoving_distance(n.arange(0,1.3,0.0001)), n.arange(0,1.3,0.0001)) 

######################
# Observations
######################
# creates a histogram of each galaxy properties, three selections

def get_hist(prop, ALL, INCLUS):
	nnInclus, bb = n.histogram(prop[INCLUS], bins=50)
	return bb, nnInclus

#######################
# MAGNETICUM SIM
#######################

top_dir = os.path.join(os.environ['OBS_REPO'], 'eRosita_lightcone_magneticum')

lc_dir = os.path.join(top_dir, 'light_cone')
#plot_dir = os.path.join(os.environ['HOME'], 'wwwDir', 'eRoMok', 'magneticum')
plot_dir = os.path.join(top_dir, 'plots')

snap_Nr, snap_z, snap_r0, snap_r1 = n.loadtxt(os.path.join(lc_dir, "wmap.geometry.dat"), unpack=True, skiprows=3)

# choose a snapshot for the clusters
id_snap = 0
id_s = str(int(snap_Nr[id_snap])).zfill(3)
print('opens cluster file', os.path.join(lc_dir, "wmap."+id_s+".cluster.fits"))
hdC = fits.open(os.path.join(lc_dir, "wmap."+id_s+".cluster.fits"))
#cl_ihal, , cl_r, cl_z_true, cl_z_obs, cl_Mvir, cl_Rvir, cl_M500, cl_R500, cl_Mstar500, cl_Mgas500, cl_T500, cl_Lx500, cl_Ysz500, cl_M200, cl_R200 = 
cl_l       = hdC[1].data['l']
cl_b       = hdC[1].data['b']
cl_z_true  = hdC[1].data['z_true']
cl_Rvir    = hdC[1].data['Rvir']

print(time.time()-t0, 'done')
#lmax = 5.
#bmax = 5.
#cl_area_sel =(cl_l<lmax)&(cl_b<bmax)

N_Rvir = 10.

rvir_min = 1800.
rvir_max = 2200.
cl_rvir_sel = (cl_Rvir > rvir_min)&(cl_Rvir < rvir_max)

cl_numbers = n.arange(len(cl_Rvir))[cl_rvir_sel]
print('N clusters', len(cl_numbers))

# get 
NN_inclus_Mstar = n.zeros((len(cl_numbers), 50))
NN_inclus_Age   = n.zeros((len(cl_numbers), 50))
NN_inclus_Z     = n.zeros((len(cl_numbers), 50))
NN_inclus_sfr   = n.zeros((len(cl_numbers), 50))
NN_inclus_g     = n.zeros((len(cl_numbers), 50))
NN_inclus_r     = n.zeros((len(cl_numbers), 50))
NN_inclus_i     = n.zeros((len(cl_numbers), 50))
NN_inclus_ug    = n.zeros((len(cl_numbers), 50))
NN_inclus_gr    = n.zeros((len(cl_numbers), 50))
NN_inclus_ri    = n.zeros((len(cl_numbers), 50))
NN_inclus_iz    = n.zeros((len(cl_numbers), 50))

# choose a cluster id : 
for id_clus in range(len(cl_numbers)):
	t1 = time.time()
	# snapshot id where the clusters are
	id_snap = 0
	id_s = str(int(snap_Nr[id_snap])).zfill(3)
	# cluster id
	cl_num = cl_numbers[id_clus]
	print(id_clus, cl_num, t1)
	# looks for objects in a cubic box around N rvir of the cluster
	print("cluster properties", cl_z_true[cl_num], cl_Rvir[cl_num])
	cl_Rvir_deg = cosmo.arcsec_per_kpc_comoving(snap_z[id_snap]).value*cl_Rvir[cl_num]/3600.
	cl_z_p_N_rvir = dC_2_z( cosmo.comoving_distance(cl_z_true[cl_num]).value + N_Rvir*cl_Rvir[cl_num]/1000.)
	cl_z_m_N_rvir = dC_2_z( cosmo.comoving_distance(cl_z_true[cl_num]).value - N_Rvir*cl_Rvir[cl_num]/1000.)
	# get the galaxies in the cluster and outside in the same snapshot slice
	print('opens galaxy file in the same snapshot as the cluster', os.path.join(lc_dir, "wmap."+id_s+".galaxies.fits"))
	#isub, l, b, rr, vmax, z_true, z_obs, Mstar, sfr, u, V,  g,  r,  i,  z,  Y,  J , H  ,K  ,L  ,M, Age ,Z = n.loadtxt(os.path.join(lc_dir, "wmap."+id_s+".galaxies.fits"), unpack=True, skiprows=1)
	hd = fits.open(os.path.join(lc_dir, "wmap."+id_s+".galaxies.fits"))
	gal_in_N_Rvir_proj = (abs(cl_l[cl_num] - hd[1].data['l'])<cl_Rvir_deg * N_Rvir ) & (abs(cl_b[cl_num] - hd[1].data['b'])<cl_Rvir_deg * N_Rvir )
	gal_in_N_Rvir_in_cluster = (gal_in_N_Rvir_proj) & (hd[1].data['z_true'] > cl_z_m_N_rvir) & (hd[1].data['z_true'] < cl_z_p_N_rvir)
	bb_Mstar, NN_inclus_Mstar[id_clus] =  get_hist(n.log10(hd[1].data['Mstar']), gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)
	bb_Age,   NN_inclus_Age[id_clus]   =  get_hist(hd[1].data['Age'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)
	bb_Z,     NN_inclus_Z[id_clus]     =  get_hist(hd[1].data['Z'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)
	bb_sfr,   NN_inclus_sfr[id_clus]   =  get_hist(n.log10(hd[1].data['sfr']+0.001), gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)
	bb_g,     NN_inclus_g[id_clus]     =  get_hist(hd[1].data['g'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)
	bb_r,     NN_inclus_r[id_clus]     =  get_hist(hd[1].data['r'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)
	bb_i,     NN_inclus_i[id_clus]     =  get_hist(hd[1].data['i'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)
	bb_ug,     NN_inclus_ug[id_clus]   =  get_hist(hd[1].data['u']-hd[1].data['g'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)
	bb_gr,     NN_inclus_gr[id_clus]   =  get_hist(hd[1].data['g']-hd[1].data['r'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)
	bb_ri,     NN_inclus_ri[id_clus]   =  get_hist(hd[1].data['r']-hd[1].data['i'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)
	bb_iz,     NN_inclus_iz[id_clus]   =  get_hist(hd[1].data['i']-hd[1].data['z'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster)


DATA = [NN_inclus_Mstar,
NN_inclus_Age  ,
NN_inclus_Z    ,
NN_inclus_sfr  ,
NN_inclus_g    ,
NN_inclus_r    ,
NN_inclus_i    ,
NN_inclus_ug    ,
NN_inclus_gr    ,
NN_inclus_ri    ,
NN_inclus_iz    ]

fb = open(os.path.join(lc_dir, "results_Nrvir_"+str(N_Rvir)+".pkl"), "wb")
pickle.dump(DATA, fb)
fb.close()

#fb = open(os.path.join(lc_dir, "results.pkl"), "r")
#DATA2 = pickle.load(fb)
#fb.close()

NN_inclus_Mstar_m = n.mean(NN_inclus_Mstar, axis=0)
NN_inclus_Age_m   = n.mean(NN_inclus_Age  , axis=0)
NN_inclus_Z_m     = n.mean(NN_inclus_Z    , axis=0)
NN_inclus_sfr_m   = n.mean(NN_inclus_sfr  , axis=0)
NN_inclus_g_m     = n.mean(NN_inclus_g    , axis=0)
NN_inclus_r_m     = n.mean(NN_inclus_r    , axis=0)
NN_inclus_i_m     = n.mean(NN_inclus_i    , axis=0)
NN_inclus_ug_m     = n.mean(NN_inclus_ug    , axis=0)
NN_inclus_gr_m     = n.mean(NN_inclus_gr    , axis=0)
NN_inclus_ri_m     = n.mean(NN_inclus_ri    , axis=0)
NN_inclus_iz_m     = n.mean(NN_inclus_iz    , axis=0)

area = 129600./n.pi/8

def plot_hist(bb, N_in, prop):
	x = (bb[1:]+bb[:-1])/2.
	p.figure(2, (8,6))
	p.axes([0.18, 0.18, 0.75, 0.75])
	p.plot(x, N_in/area, label='members', ls='solid', lw=3)
	p.xlabel(prop)
	p.yscale('log')
	p.ylabel('Number in 10 Rvir / deg2')
	p.title(r'clusters with $1.8<rvir/Mpc<2.2$ at $z\sim0.25$')
	p.legend(frameon=False)
	p.savefig(os.path.join(plot_dir, 'hist_'+str(N_Rvir)+'_'+prop+'.png'))
	p.clf()

plot_hist(bb_ug, NN_inclus_ug_m, "U-G ABS MAG")
plot_hist(bb_gr, NN_inclus_gr_m, "G-R ABS MAG")
plot_hist(bb_ri, NN_inclus_ri_m, "R-I ABS MAG")
plot_hist(bb_iz, NN_inclus_iz_m, "I-Z ABS MAG")
plot_hist(bb_g, NN_inclus_g_m,   "G ABS MAG")
plot_hist(bb_r, NN_inclus_r_m,   "R ABS MAG")
plot_hist(bb_i, NN_inclus_i_m,   "I ABS MAG")
plot_hist(bb_Mstar, NN_inclus_Mstar_m, 'Stellar mass')
plot_hist(bb_Age, NN_inclus_Age_m   , 'Stellar Age')
plot_hist(bb_Z, NN_inclus_Z_m     ,'Stellar metallicity')
plot_hist(bb_sfr, NN_inclus_sfr_m   , 'Star formation rate')




