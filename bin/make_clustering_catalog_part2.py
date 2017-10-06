#!/usr/bin/env python

import astropy.io.fits as fits
import os
import sys
from os.path import join

import numpy as n
from scipy.interpolate import interp1d
import scipy.spatial.ckdtree as t

path_2_plate_file = join( os.environ['SDSS_DATA_DIR'], 'dr14','plates-dr14.fits')
plates = fits.open(path_2_plate_file)[1].data
target_dir = join( os.environ['SDSS_DATA_DIR'], "dr14", "spiders", "target")

path_2_target_AGN_RASS_eboss  = os.path.join(target_dir, "SPIDERS_RASS_AGN_targets_for_EBOSS_TDSS_SPIDERS.fits.gz")
path_2_target_AGN_RASS_sequels  = os.path.join(target_dir, "SPIDERS_RASS_AGN_targets_for_SEQUELS.fits.gz")
path_2_target_AGN_XMMSL_eboss  = os.path.join(target_dir, "SPIDERS_XMMSL_AGN_targets_for_EBOSS_TDSS_SPIDERS.fits.gz")



r_cp = 92./3600. # degrees, center post mask

# thresholds to construct the catalogs
w_min = 0.2 # 30% sampling rate
n_min = 1. # minimum of 1 objects in a polygon
z_min = 0.
z_max = 4.

def create_clustering_catalogs(spec_file, rd_file, completeness_per_polyid_file, w_min=w_min, n_min=n_min):
	print os.path.isfile(spec_file), spec_file
	print os.path.isfile(rd_file), rd_file
	print os.path.isfile(completeness_per_polyid_file), completeness_per_polyid_file
	output_phot_cat = spec_file[:-5]+"_clustering.fits"
	ra_tile, dec_tile = plates['RACEN'], plates['DECCEN']
	print " first the data "
	hd_obs     = fits.open(spec_file)[1].data
	treeD = t.cKDTree( n.transpose([hd_obs['ra'], hd_obs['dec']]) )
	treeT = t.cKDTree( n.transpose([ra_tile, dec_tile]) )
	id_cp = treeT.query_ball_tree(treeD, r_cp)
	in_cp_id_list = n.array(list(set(n.hstack((id_cp)))),dtype=n.int32)
	print "in cp mask", len(in_cp_id_list)
	in_cp_id_list.sort()
	cp_mask = n.ones(len(hd_obs['ra']))
	cp_mask[in_cp_id_list] = n.zeros_like(cp_mask[in_cp_id_list])
	sample = (hd_obs['nin']>n_min)&(hd_obs['TSR']*hd_obs['SSR']>w_min)&(hd_obs['Z_ERR']>0)&(hd_obs['Z']>z_min)&(hd_obs['Z']<z_max)
	w_col = fits.Column(name="w",format='D', array=hd_obs['TSR'][sample]*hd_obs['SSR'][sample] )
	ra_col = fits.Column(name="ra",format='D', array=hd_obs['ra'][sample])
	dec_col = fits.Column(name="dec",format='D', array=hd_obs['dec'][sample])
	z_col = fits.Column(name="z",format='D', array=hd_obs['Z'][sample])
	print "z range data: ", n.min(z_col.array), n.max(z_col.array)
	w_err_col = fits.Column(name="w_err",format='D', array=hd_obs['TSR'][sample]*hd_obs['SSR'][sample]*hd_obs['nin'][sample]**(-0.5) )
	z_err_col = fits.Column(name="z_err",format='D', array=hd_obs['Z_ERR'][sample])
	coldefs = fits.ColDefs([ra_col, dec_col, z_col, w_col, z_err_col, w_err_col])
	hdu = fits.BinTableHDU.from_columns(coldefs)
	if os.path.isfile(output_phot_cat):
		os.remove(output_phot_cat)

	hdu.writeto(output_phot_cat)
	print " now the randoms "
	all_polys, nin_poly, tsr_poly, ssr_poly = n.loadtxt(completeness_per_polyid_file, unpack=True)
	poly_ok = (nin_poly > n_min ) & (tsr_poly * ssr_poly> w_min)
	polys = all_polys[poly_ok]
	rds_outfile = rd_file[:-5]+"_clustering.fits"
	hd_rds     = fits.open(rd_file)[1].data
	treeR = t.cKDTree( n.transpose([hd_rds['ra'], hd_rds['dec']]) )
	treeT = t.cKDTree( n.transpose([ra_tile, dec_tile]) )
	id_cp = treeT.query_ball_tree(treeR, r_cp)
	in_cp_id_list = n.array(list(set(n.hstack((id_cp)))),dtype=n.int32)
	print "randoms in cp mask", len(in_cp_id_list)
	in_cp_id_list.sort()
	cp_mask = n.ones(len(hd_rds['ra']))
	cp_mask[in_cp_id_list] = n.zeros_like(cp_mask[in_cp_id_list])
	unmask = (cp_mask==1)
	hd_rd_inter = n.hstack((n.array([hd_rds[(unmask)&(hd_rds['tiling_polyid']==poly)] for poly in polys])))
	nR = len(hd_rd_inter)
	dz=0.05
	zs=n.arange(0.,z_col.array.max()+dz,dz)
	nn,bb = n.histogram(z_col.array, bins=zs, weights=1./w_col.array)
	nz=interp1d((zs[1:]+zs[:-1])/2.,nn)
	rdsz=[]
	for i in range(1,len(zs)-1,1):
		inter=n.random.uniform(low=zs[i]-dz/2., high=zs[i]+dz/2., size=int( 1000* nz( zs[i] )))
		rdsz.append(inter)

	rds=n.hstack((rdsz))
	n.random.shuffle(rds)
	RR=rds[:nR]#-dz/2.

	ra_rd_col = fits.Column(name="ra",format='D', array=hd_rd_inter['RA'])
	dec_rd_col = fits.Column(name="dec",format='D', array=hd_rd_inter['DEC'])
	z_rd_col = fits.Column(name="z",format='D', array=RR)

	coldefs = fits.ColDefs([ra_rd_col,  dec_rd_col, z_rd_col])
	hdu = fits.BinTableHDU.from_columns(coldefs)

	if os.path.isfile(rds_outfile):
		os.remove(rds_outfile)

	hdu.writeto(rds_outfile)



spec_file = path_2_target_AGN_XMMSL_eboss[:-8]+"_merged_TSR.fits"
rd_file = path_2_target_AGN_XMMSL_eboss[:-8]+"_merged_random.fits"
completeness_per_polyid_file = path_2_target_AGN_XMMSL_eboss[:-8]+"_merged_completeness.txt"
create_clustering_catalogs(spec_file, rd_file, completeness_per_polyid_file, w_min=w_min, n_min=n_min)

spec_file = path_2_target_AGN_RASS_eboss[:-8]+"_merged_TSR.fits"
rd_file = path_2_target_AGN_RASS_eboss[:-8]+"_merged_random.fits"
completeness_per_polyid_file = path_2_target_AGN_RASS_eboss[:-8]+"_merged_completeness.txt"
create_clustering_catalogs(spec_file, rd_file, completeness_per_polyid_file, w_min=w_min, n_min=n_min)

spec_file = path_2_target_AGN_RASS_sequels[:-8]+"_merged_TSR.fits"
rd_file = path_2_target_AGN_RASS_sequels[:-8]+"_merged_random.fits"
completeness_per_polyid_file = path_2_target_AGN_RASS_sequels[:-8]+"_merged_completeness.txt"
create_clustering_catalogs(spec_file, rd_file, completeness_per_polyid_file, w_min=w_min, n_min=n_min)

