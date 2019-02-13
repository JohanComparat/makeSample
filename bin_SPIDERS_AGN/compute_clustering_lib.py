import numpy as n
import lc_setup
import astropy.io.fits as fits
import h5py    # HDF5 support
import os
import glob
import numpy as n
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)

theory_dir = '/data17s/darksim/theory/camb'
		

# catalogs :
# /data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/2RXS/2RXS_AllWISE_catalog_paper_2017May26.fits.gz
# /data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/XMM/XMMSL2_AllWISE_catalog_paper_2017JUN09.fits.gz
# /data44s/eroAGN_WG_DATA/DATA/randoms/randoms.fits

# sub catalogs
# - full cat
# - cut abs(gal_lat)>20 degrees
# - cut on X-ray depth map
# - cut on wise depth map

# write ascii catalogs + param files here :
# 

out_dir='/data36s/comparat/AGN_clustering/angular_clustering/'

def get_ra_dec_z_data(data, selection, out_name ):
	out_name = os.path.join(out_dir , out_name + '.data')
	n_gal = len( data['RA'][selection] )
	n.savetxt(out_name, n.transpose([data['RA'][selection], data['DEC'][selection], 0.7*n.ones_like(data['RA'][selection]), n.ones_like(data['RA'][selection]) ])  )
	print(out_name, n_gal)
	return n_gal

def create_randoms(Lname, tracer_name, zmin, zmax):
	out_dir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_'+Lname+'/'
	out_name = os.path.join(out_dir , Lname + '_' + tracer_name + '_' + str(zmin) + '_' + str(zmax) + '.random')
	print('loads data')
	raD, decD, zD, wD = n.loadtxt(os.path.join(out_dir , Lname + '_' + tracer_name + '_' + str(zmin) + '_' + str(zmax) +'.data'), unpack=True)
	print('opens randoms')
	raR, decR = n.loadtxt(out_dir + 'random-ra-dec.txt', unpack=True)
	# the random
	N_data = len(raD) 
	N_rds = 10*N_data # len(raR) 
	print("D,R=",N_data, N_rds)
	dz=0.025
	zs=n.arange(zmin, zmax + 2*dz, dz)
	nn,bb = n.histogram(zD, bins=zs)#, weights=1./w_col.array)
	#nz=interp1d((zs[1:]+zs[:-1])/2.,nn)
	rdsz=[]
	for jj, (zlow, zhigh) in enumerate(zip(zs[:-1], zs[1:])):
		#print(zlow, zhigh)
		inter=n.random.uniform(low=zlow, high=zhigh, size=int( 1000* nn[jj] ))
		rdsz.append(inter)

	rds=n.hstack((rdsz))
	n.random.shuffle(rds)
	RR=rds[:N_rds]#-dz/2.
	print("RR=",len(rds), len(RR))
	n.random.shuffle(raR)
	n.random.shuffle(decR)
	n.savetxt(out_name, n.transpose([raR[:N_rds], decR[:N_rds], RR, n.ones_like(RR) ]))

def compute_clustering(Lname,  tracer_name, zmin, zmax):
	out_dir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_'+Lname+'/'
	out_name = os.path.join(out_dir , Lname + '_' + tracer_name + '_' + str(zmin) + '_' + str(zmax) )
	f=open(out_name+'.ini','w')
	f.write('data_filename= '+out_name+'.data \n')    
	f.write('random_filename= '+out_name+'.random \n')    
	f.write('input_format= 2 \n')
	f.write('output_filename= '+out_name+'.2pcf \n')
	f.write('corr_type= monopole \n')
	f.write('omega_M= 0.307 \n')
	f.write('omega_L= 0.693 \n')
	f.write('w= -1 \n')
	f.write('log_bin= 0 \n')
	f.write('dim1_max= 160. \n')
	#f.write('dim1_min_logbin= 0.01 \n')
	f.write('dim1_nbin= 40 \n')
	f.write('dim2_max= 160. \n')
	f.write('dim2_nbin= 40 \n')
	f.write('dim3_min= 0.00 \n')
	f.write('dim3_max= 3. \n')
	f.write('dim3_nbin= 1 \n')
	f.close()
	os.system("~/darksim/software/CUTE/CUTE/CUTE "+out_name+'.ini')


def plot_results(Lname,  tracer_name, zmin, zmax):
	out_dir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_'+Lname+'/'
	out_name = os.path.join(out_dir , Lname + '_' + tracer_name + '_' + str(zmin) + '_' + str(zmax) )
	setup = lc_setup.cones[Lname]
	area = setup.area
	dz=0.05
	z_bins = n.arange(zmin, zmax+dz, dz)
	xi_file = 'xi_0.3.dat'
	dec_max = 3.	

	fig = p.figure(2, (10,8))
	# wedge plot
	fig.add_subplot(2,1,1)

	ra, dec, z, w = n.loadtxt(out_name + '.data', unpack=True)
	ok = (abs(dec)<dec_max)
	RR_all = cosmoMD.comoving_distance(z[ok]).value
	x_all = RR_all*n.cos(ra[ok]* n.pi/180.) 
	y_all = RR_all*n.sin(dec[ok]* n.pi/180.)
	p.plot(x_all, y_all, marker='.', color='b', alpha=0.1, rasterized = True, ls='None')
	p.xlabel('x [Mpc]')
	p.ylabel('y [Mpc]')
	p.title(os.path.basename(out_name))

	fig.add_subplot(2,2,3)
	# nz
	p.hist(z, bins = z_bins, weights=n.ones_like(z)/area)
	p.xlabel('redshift')
	p.ylabel('N/deg2')
	fig.add_subplot(2,2,4)
	# clustering
	rs,xs = n.loadtxt(os.path.join(theory_dir, xi_file), unpack=True)
	p.plot(rs[::3], rs[::3]*rs[::3]*xs[::3]*2 )#, label=label)
	s420 = n.loadtxt(out_name+'.2pcf' , unpack=True)
	p.plot(s420 [0], s420 [0]**2*s420 [1])#, label=label)
	p.xlabel('s [Mpc/h]')
	p.ylabel(r'$s^2.\xi(s)$')
	#p.xscale('log')
	#p.yscale('log')
	p.ylim((-20,200))
	p.xlim((0.5,165))

	p.savefig(out_name+".png")
	p.clf()








