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

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as p

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

"""
path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Finoguenov_cosmos_2007_clusters.data')
x_data, y_data, y_min, y_max = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.fill_between(n.log10(x_data), y1 = n.log10(y_max), y2=n.log10(y_min), color='b' , rasterized = True, alpha=0.5)
p.plot(n.log10(x_data), n.log10(y_data), color='b', label = 'COSMOS Finoguenov 2007' )

path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Finoguenov_ecdfs_2015_clusters.data')
x_data, y_data, y_min, y_max = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.fill_between(n.log10(x_data), y1 = n.log10(y_max), y2=n.log10(y_min) , rasterized = True, alpha=0.5, color='g' )
p.plot(n.log10(x_data), n.log10(y_data), color='g',  label = 'ECDFS Finoguenov 2015' )
"""


# creates a histogram of each galaxy properties, three selections

def get_hist(prop, ALL, INCLUS, CONTAMINATION):
	nnInclus, bb = n.histogram(prop[INCLUS], bins=50)
	nnContamination = n.histogram(prop[CONTAMINATION], bins=bb)[0]
	return bb, nnInclus, nnContamination


def get_hist_1(prop, CONTAMINATION, bins):
	nnContamination = n.histogram(prop[CONTAMINATION], bins=bins)[0]
	return nnContamination


# get the galaxies in the cluster and outside in a different snapshot slice

def get_all_hist(lc_dir, id_snap, cl_l, cl_b, cl_num, cl_Rvir_deg, N_Rvir, bb_Mstar, bb_Age, bb_Z, bb_sfr, bb_g, bb_r, bb_i):
	id_s = str(int(snap_Nr[id_snap])).zfill(3)
	hdu = fits.open(os.path.join(lc_dir, "wmap."+id_s+".galaxies.fits"))
	#isub, l, b, rr, vmax, z_true, z_obs, Mstar, sfr, u, V,  g,  r,  i,  z,  Y,  J , H  ,K  ,L  ,M, Age ,Z = 
	gal_in_N_Rvir_proj = (abs(cl_l[cl_num] - hdu[1].data['l'])<cl_Rvir_deg * N_Rvir ) & (abs(cl_b[cl_num] - hdu[1].data['b'])<cl_Rvir_deg * N_Rvir )
	NN_contamination_Mstar =  get_hist_1(n.log10(hdu[1].data['Mstar']),  gal_in_N_Rvir_proj, bb_Mstar)
	NN_contamination_Age   =  get_hist_1(hdu[1].data['Age'], gal_in_N_Rvir_proj, bb_Age)
	NN_contamination_Z     =  get_hist_1(hdu[1].data['Z'], gal_in_N_Rvir_proj, bb_Z)
	NN_contamination_sfr   =  get_hist_1(n.log10(hdu[1].data['sfr']+0.001), gal_in_N_Rvir_proj, bb_sfr)
	NN_contamination_g     =  get_hist_1(hdu[1].data['g'], gal_in_N_Rvir_proj, bb_g)
	NN_contamination_r     =  get_hist_1(hdu[1].data['r'], gal_in_N_Rvir_proj, bb_r)
	NN_contamination_i     =  get_hist_1(hdu[1].data['i'], gal_in_N_Rvir_proj, bb_i)
	return NN_contamination_Mstar, NN_contamination_Age, NN_contamination_Z, NN_contamination_sfr, NN_contamination_g, NN_contamination_r, NN_contamination_i  

#######################
# MAGNETICUM SIM
#######################

top_dir = os.path.join(os.environ['OBS_REPO'], 'eRosita_lightcone_magneticum')

lc_dir = os.path.join(top_dir, 'light_cone')
plot_dir = os.path.join(os.environ['HOME'], 'wwwDir', 'eRoMok', 'magneticum')

snap_Nr, snap_z, snap_r0, snap_r1 = n.loadtxt(os.path.join(lc_dir, "wmap.geometry.dat"), unpack=True, skiprows=3)

# choose a snapshot for the clusters
id_snap = -1
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

N_Rvir = 5.

rvir_min = 1800.
rvir_max = 2200.
cl_rvir_sel = (cl_Rvir > rvir_min)&(cl_Rvir < rvir_max)

cl_numbers = n.arange(len(cl_Rvir))[cl_rvir_sel]
print('N clusters', len(cl_numbers))


NN_contamination_Mstar = n.zeros((len(cl_numbers), len(snap_z), 50))
NN_contamination_Age   = n.zeros((len(cl_numbers), len(snap_z), 50))
NN_contamination_Z     = n.zeros((len(cl_numbers), len(snap_z), 50))
NN_contamination_sfr   = n.zeros((len(cl_numbers), len(snap_z), 50))
NN_contamination_g     = n.zeros((len(cl_numbers), len(snap_z), 50))
NN_contamination_r     = n.zeros((len(cl_numbers), len(snap_z), 50))
NN_contamination_i     = n.zeros((len(cl_numbers), len(snap_z), 50))

NN_inclus_Mstar = n.zeros((len(cl_numbers), 50))
NN_inclus_Age   = n.zeros((len(cl_numbers), 50))
NN_inclus_Z     = n.zeros((len(cl_numbers), 50))
NN_inclus_sfr   = n.zeros((len(cl_numbers), 50))
NN_inclus_g     = n.zeros((len(cl_numbers), 50))
NN_inclus_r     = n.zeros((len(cl_numbers), 50))
NN_inclus_i     = n.zeros((len(cl_numbers), 50))


# choose a cluster id : 
for id_clus in range(len(cl_numbers)):
	t1 = time.time()
	# snapshot id where the clusters are
	id_snap = -1
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
	gal_in_N_Rvir_contamination = (gal_in_N_Rvir_proj) & (gal_in_N_Rvir_in_cluster==False)
	#print("in N Rvir (all, in cluster, contamination):", len(z_true[gal_in_N_Rvir_proj]),len(z_true[gal_in_N_Rvir_in_cluster]), len(z_true[gal_in_N_Rvir_contamination]))
	bb_Mstar, NN_inclus_Mstar[id_clus], NN_contamination_Mstar[id_clus][id_snap]  =  get_hist(n.log10(hd[1].data['Mstar']), gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster, gal_in_N_Rvir_contamination)
	bb_Age,   NN_inclus_Age[id_clus],   NN_contamination_Age[id_clus][id_snap]    =  get_hist(hd[1].data['Age'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster, gal_in_N_Rvir_contamination)
	bb_Z,     NN_inclus_Z[id_clus],     NN_contamination_Z[id_clus][id_snap]      =  get_hist(hd[1].data['Z'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster, gal_in_N_Rvir_contamination)
	bb_sfr,   NN_inclus_sfr[id_clus],   NN_contamination_sfr[id_clus][id_snap]    =  get_hist(n.log10(hd[1].data['sfr']+0.001), gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster, gal_in_N_Rvir_contamination)
	bb_g,     NN_inclus_g[id_clus],     NN_contamination_g[id_clus][id_snap]      =  get_hist(hd[1].data['g'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster, gal_in_N_Rvir_contamination)
	bb_r,     NN_inclus_r[id_clus],     NN_contamination_r[id_clus][id_snap]      =  get_hist(hd[1].data['r'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster, gal_in_N_Rvir_contamination)
	bb_i,     NN_inclus_i[id_clus],     NN_contamination_i[id_clus][id_snap]      =  get_hist(hd[1].data['i'], gal_in_N_Rvir_proj, gal_in_N_Rvir_in_cluster, gal_in_N_Rvir_contamination)
		#print(time.time()-t0, 'histograms done')
	for id_snap in n.arange(1,len(snap_Nr),1)[::-1]:
		NN_contamination_Mstar[id_clus][id_snap], NN_contamination_Age[id_clus][id_snap], NN_contamination_Z[id_clus][id_snap], NN_contamination_sfr[id_clus][id_snap], NN_contamination_g[id_clus][id_snap], NN_contamination_r[id_clus][id_snap], NN_contamination_i[id_clus][id_snap] = get_all_hist(lc_dir, id_snap, cl_l, cl_b, cl_num, cl_Rvir_deg, N_Rvir, bb_Mstar, bb_Age, bb_Z, bb_sfr, bb_g, bb_r, bb_i)
		print('snap', id_snap, time.time()-t1)


DATA = n.array([NN_contamination_Mstar ,
NN_contamination_Age   ,
NN_contamination_Z     ,
NN_contamination_sfr   ,
NN_contamination_g     ,
NN_contamination_r     ,
NN_contamination_i     ,
NN_inclus_Mstar,
NN_inclus_Age  ,
NN_inclus_Z    ,
NN_inclus_sfr  ,
NN_inclus_g    ,
NN_inclus_r    ,
NN_inclus_i    ])

pickle.dump(DATA, open(os.path.join(lc_dir, "results.pkl"), "w"))

#def plot_hist(bb, N_in, N_cont, prop):
	#x = (bb[1:]+bb[:-1])/2.
	#p.figure(0, (5,5))
	#p.axes([0.18, 0.18, 0.75, 0.75])
	#p.plot(x, N_in, label='members', ls='solid', lw=3)
	#for ii, el in enumerate(N_cont):
		#p.plot(x, el, label='z='+str(snap_z[ii]), ls='dashed')

	#p.xlabel(prop)
	#p.ylabel('Number in 5 Rvir')
	#p.savefig(os.path.savefig(plot_dir, 'hist_'+prop+'.png'))
	#p.clf()

#plot_hist(bb_g, NN_inclus_g, NN_contamination_g, "G ABS MAG")
#plot_hist(bb_r, NN_inclus_r, NN_contamination_r, "R ABS MAG")
#plot_hist(bb_i, NN_inclus_i, NN_contamination_i, "I ABS MAG")