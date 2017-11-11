# python 2 to have healpy
from os.path import join
import os
import glob
import time

#from _pickle import cPickle
import cPickle
#import fileinput
import astropy.io.fits as fits

import numpy as n

import healpy 

import astropy.cosmology as co
import astropy.units as u
cosmo = co.WMAP7

flux_bins = n.arange(-18, -8., 0.2)

#######################
# MAGNETICUM SIM
#######################

top_dir = os.path.join(os.environ['OBS_REPO'], 'eRosita_lightcone_magneticum')

lc_dir = os.path.join(top_dir, 'light_cone')
plot_dir = os.path.join(os.environ['HOME'], 'wwwDir', 'eRoMok', 'magneticum')

snap_Nr, snap_z, snap_r0, snap_r1 = n.loadtxt(os.path.join(lc_dir, "wmap.geometry.dat"), unpack=True, skiprows=3)

# choose a snapshot for the clusters
id_snap = 0
def get_fluxes(id_snap):
	id_s = str(int(snap_Nr[id_snap])).zfill(3)
	print('opens cluster file', os.path.join(lc_dir, "wmap."+id_s+".cluster.fits"))
	hdC = fits.open(os.path.join(lc_dir, "wmap."+id_s+".cluster.fits"))
	#cl_ihal, , cl_r, cl_z_true, cl_z_obs, cl_Mvir, cl_Rvir, cl_M500, cl_R500, cl_Mstar500, cl_Mgas500, cl_T500, cl_Lx500, cl_Ysz500, cl_M200, cl_R200 = 
	#cl_l       = hdC[1].data['l']
	#cl_b       = hdC[1].data['b']
	#cl_Rvir    = hdC[1].data['Rvir']
	cl_z_true  = hdC[1].data['z_true']
	cl_LX_500  = hdC[1].data['Lx500'] * 1e44
	dl_cm = cosmo.luminosity_distance(cl_z_true).to(u.cm)
	cl_fl_500 = cl_LX_500/(4*n.pi*dl_cm**2.)
	return cl_fl_500

fl_500 = []
for id_snap in n.arange(len(snap_Nr)):
	fl_500.append(get_fluxes(id_snap))

magneticum_fl500 = n.hstack(fl_500)
magneticum_area = 129600./n.pi/8.

out_magneticum = n.histogram(n.log10(magneticum_fl500), bins = flux_bins)[0]
c_out_magneticum = n.array([n.sum(out_magneticum[ii:]) for ii in range(len(out_magneticum)) ]) / magneticum_area

#ras = n.array([healpy.pix2ang(NSIDE,pix_id)[1]*180./n.pi for pix_id in pix_ids ])
#decs = n.array([ (healpy.pix2ang(NSIDE,pix_id)[0]-n.pi/2.)*180./n.pi for pix_id in pix_ids ])
#n.savetxt("px_ra_Dec.txt", n.transpose([pix_ids, ras, decs]))
NSIDE = 16
pix_ids_16 = n.arange(healpy.nside2npix(NSIDE))
area_per_pixel_16 = 129600./n.pi/healpy.nside2npix(NSIDE)
print("pixel area considered=",area_per_pixel_16,"deg2")

NSIDE = 32
pix_ids_32 = n.arange(healpy.nside2npix(NSIDE))
area_per_pixel_32 = 129600./n.pi/healpy.nside2npix(NSIDE)
print("pixel area considered=",area_per_pixel_32,"deg2")


path_2_light_cone = os.path.join(os.environ["MD10"], 'light-cone', 'MDPL2_FluxProj000_ClustersCombinedModel_bias0.6_with_header.fits')
hd = fits.open(path_2_light_cone)
log_f_05_20 = n.log10(hd[1].data['F_05_20'])

HEALPIX_16 = healpy.ang2pix(16,  hd[1].data['galactic_latitude_deg']*n.pi/180. , hd[1].data['galactic_longitude_deg']*n.pi/180. ) 
HEALPIX_32 = healpy.ang2pix(32,  hd[1].data['galactic_latitude_deg']*n.pi/180. , hd[1].data['galactic_longitude_deg']*n.pi/180. ) 

out = n.histogram(log_f_05_20, bins = flux_bins)

per_pixel_out_16 = n.array([ n.histogram(log_f_05_20[HEALPIX_16==hp_i], bins = n.arange(-18, -8., 0.2))[0] for hp_i in pix_ids_16[:1000] ]) 
per_pixel_out_c_16 = n.array([[ n.sum(el[ii:]) for ii in range(len(el)) ] for el in per_pixel_out_16 ])
frac_err_13deg2 = n.std(per_pixel_out_c_16, axis=0)/n.mean(per_pixel_out_c_16, axis=0)

per_pixel_out_32 = n.array([ n.histogram(log_f_05_20[HEALPIX_32==hp_i], bins = n.arange(-18, -8., 0.2))[0] for hp_i in pix_ids_32[:1000] ]) 
per_pixel_out_c_32 = n.array([[ n.sum(el[ii:]) for ii in range(len(el)) ] for el in per_pixel_out_32 ])
frac_err_3deg2 = n.std(per_pixel_out_c, axis=0)/n.mean(per_pixel_out_c_32, axis=0)

# cumulative number density per square degrees
x_out = 0.5*(out[1][1:] + out[1][:-1])
c_out = n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ])*n.pi/129600.


#path_2_light_cone = os.path.join(os.environ["MD10"], 'light-cone', 'MDPL2_FluxProj000_ClustersCombinedModel_bias0.6_with_header.fits')
#hd = fits.open(path_2_light_cone)
#log_f_05_20 = n.log10(hd[1].data['F_05_20'])
#out = n.histogram(log_f_05_20, bins = n.arange(-18, -8., 0.2))
## cumulative number density per square degrees
#x_out = 0.5*(out[1][1:] + out[1][:-1])
#c_out = n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ])*n.pi/129600.

#path_2_light_cone = os.path.join(os.environ["MD10"], 'light-cone', 'MDPL2_FluxProj000_ClustersCombinedModel_with_header.fits')
#hd = fits.open(path_2_light_cone)
#log_f_05_20 = n.log10(hd[1].data['F_05_20'])
#out_0 = n.histogram(log_f_05_20, bins = n.arange(-18, -8., 0.2))
## cumulative number density per square degrees
#x_out_0 = 0.5*(out_0[1][1:] + out_0[1][:-1])
#c_out_0 = n.array([n.sum(out_0[0][ii:]) for ii in range(len(out_0[0])) ])*n.pi/129600.


path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Finoguenov_07_15_clusters.data')
x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

plotDir = os.path.join(os.environ['HOME'], 'wwwDir', "eRoMok", "logNlogS")

p.figure(1, (6,6))
p.plot(x_out, n.log10(c_out), 'k', lw=2, rasterized = True, label = 'Planck mock v0.6' )
p.plot(x_out, n.log10(c_out*(1-frac_err_13deg2)), 'k--', lw=1, rasterized = True, label = 'v0.6, 13.3deg2 scatter' )
p.plot(x_out, n.log10(c_out*(1+frac_err_13deg2)), 'k--', lw=1, rasterized = True)
p.plot(x_out, n.log10(c_out*(1-frac_err_3deg2)), 'r--', lw=1, rasterized = True, label = 'v0.6, 3.5deg2 scatter' )
p.plot(x_out, n.log10(c_out*(1+frac_err_3deg2)), 'r--', lw=1, rasterized = True)
p.plot(x_out, n.log10(c_out_magneticum), 'm', rasterized = True, label = 'Magneticum' )

path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Finoguenov_cosmos_2007_clusters.data')
x_data, y_data, y_min, y_max = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.fill_between(n.log10(x_data), y1 = n.log10(y_max), y2=n.log10(y_min), color='b' , rasterized = True, alpha=0.5)
p.plot(n.log10(x_data), n.log10(y_data), color='b', label = 'COSMOS Finoguenov 2007' )
path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Finoguenov_ecdfs_2015_clusters.data')
x_data, y_data, y_min, y_max = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.fill_between(n.log10(x_data), y1 = n.log10(y_max), y2=n.log10(y_min) , rasterized = True, alpha=0.5, color='g' )
p.plot(n.log10(x_data), n.log10(y_data), color='g',  label = 'ECDFS Finoguenov 2015' )

p.axhline(7, ls='dashed')
p.xlabel('log(F[0.5-2 keV])')
p.ylabel('log(>F) [/deg2]')
p.legend(frameon=False, loc=0)
#p.yscale('log')
p.xlim((-16, -12))
p.ylim((-2, 3.1))
p.title('cluster mocks')
p.grid()
p.savefig(os.path.join(plotDir, "logN_logS_clusters.jpg"))
p.clf()


p.figure(1, (6,6))
p.plot(10**x_out, c_out, 'k', lw=2, rasterized = True, label = 'Planck mock v0.6' )
p.plot(10**x_out, c_out*(1-frac_err_13deg2), 'k--', lw=1, rasterized = True, label = 'v0.6, 13.3deg2 scatter' )
p.plot(10**x_out, c_out*(1+frac_err_13deg2), 'k--', lw=1, rasterized = True)
p.plot(10**x_out, c_out*(1-frac_err_3deg2), 'r--', lw=1, rasterized = True, label = 'v0.6, 3.5deg2 scatter' )
p.plot(10**x_out, c_out*(1+frac_err_3deg2), 'r--', lw=1, rasterized = True)
p.plot(x_out, n.log10(c_out_magneticum), 'm', rasterized = True, label = 'Magneticum' )

path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Finoguenov_cosmos_2007_clusters.data')
x_data, y_data, y_min, y_max = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.fill_between(x_data, y1 = y_max, y2=y_min, color='b' , rasterized = True, alpha=0.5)
p.plot(x_data, y_data, color='b', label = 'COSMOS Finoguenov 2007' )

path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Finoguenov_ecdfs_2015_clusters.data')
x_data, y_data, y_min, y_max = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.fill_between(x_data, y1 = y_max, y2=y_min , rasterized = True, alpha=0.5, color='g' )
p.plot(x_data, y_data, color='g',  label = 'ECDFS Finoguenov 2015' )

p.axhline(7, ls='dashed')
p.xlabel('F[0.5-2 keV]')
p.ylabel('>F [/deg2]')
p.legend(frameon=False, loc=0)
p.yscale('log')
p.xscale('log')
p.xlim((1e-16, 1e-12))
p.ylim((1e-2, 2e3))
p.title('cluster mocks')
p.grid()
p.savefig(os.path.join(plotDir, "logN_logS_clusters_loglog.jpg"))
p.clf()















