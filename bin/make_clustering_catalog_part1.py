#!/usr/bin/env python

import astropy.io.fits as fits
import os
import sys
from os.path import join
import numpy as n
import scipy.spatial.ckdtree as t
import pymangle as mangle

#import matplotlib.pyplot as p

spiders_dir = join( os.environ['SDSS_DATA_DIR'], "dr14", "spiders")

target_dir = join( os.environ['SDSS_DATA_DIR'], "dr14", "spiders", "target")
spectro_dir = join( os.environ['SDSS_DATA_DIR'], "dr14", "spiders", "spectroscopy")
analysis_dir = join( os.environ['SDSS_DATA_DIR'], "dr14", "spiders", "analysis")

path_2_VAC_spiders_2RXS_DR14 = os.path.join(analysis_dir, "VAC_spiders_2RXS_DR14.fits")
path_2_VAC_spiders_XMMSL_DR14 = os.path.join(analysis_dir, "VAC_spiders_XMMSL_DR14.fits")
path_2_RASS_sources_with_DR12_spectra  = os.path.join(spectro_dir, "SPIDERS_RASS_sources_with_DR12_spectra.fits.gz")
path_2_XMMSL_sources_with_DR12_spectra = os.path.join(spectro_dir, "SPIDERS_XMMSL_sources_with_DR12_spectra.fits.gz")

#path_2_target_AGN_RASS_eboss  = os.path.join(target_dir, "spiderstargetAGN-SPIDERS_RASS_AGN-v2.1.fits")
#path_2_target_AGN_RASS_sequels  = os.path.join(target_dir, "spiderstargetSequelsAGN-SPIDERS_RASS_AGN-v1.1.fits")
#path_2_target_AGN_XMMSL_eboss  = os.path.join(target_dir, "spiderstargetAGN-SPIDERS_XMMSL_AGN-v3.1.fits")
#alternative catalogs 
path_2_target_AGN_RASS_eboss  = os.path.join(target_dir, "SPIDERS_RASS_AGN_targets_for_EBOSS_TDSS_SPIDERS.fits.gz")
path_2_target_AGN_RASS_sequels  = os.path.join(target_dir, "SPIDERS_RASS_AGN_targets_for_SEQUELS.fits.gz")
path_2_target_AGN_XMMSL_eboss  = os.path.join(target_dir, "SPIDERS_XMMSL_AGN_targets_for_EBOSS_TDSS_SPIDERS.fits.gz")

path_2_specObjAll = os.path.join( os.environ['SDSS_DATA_DIR'], "dr14", "specObj-BOSS-dr14.fits")

# loads the DR14 mask
#path_2_mask_fits = os.path.join(os.environ['SDSS_DATA_DIR'], 'LSS', 'mask-lrg-N-eboss_v1.8_IRt.fits')
path_2_mask_ply = os.path.join(os.environ['SDSS_DATA_DIR'], 'LSS', 'LRG', 'v1.8', 'mask-lrg-N-eboss_v1.8_IRt.ply')
mng         = mangle.Mangle(path_2_mask_ply)


def construct_summary_files_data(spec_file, phot_file, geometry_file, output_phot_cat, rmax = 0.5/3600. ):
	# spectroscopy
	print os.path.isfile(spec_file), spec_file
	print os.path.isfile(phot_file), phot_file
	print os.path.isfile(geometry_file), geometry_file
	tiling_mask = mangle.Mangle(geometry_file)

	hd_p = fits.open(phot_file)[1].data
	hd_s_i = fits.open(spec_file)[1].data
	hd_s = hd_s_i[hd_s_i['ZWARNING']==0]

	treeP = t.cKDTree( n.transpose([hd_p['ra'], hd_p['dec']]) )
	treeS = t.cKDTree( n.transpose([hd_s['plug_ra'], hd_s['plug_dec']]) )

	Ntarg = len(hd_p)
	Nspec = len(hd_s)

	print "N targeted in chunk = ", Ntarg

	# Adding columns to the complete photometry file
	polyid_val_p      = tiling_mask.polyid(hd_p['ra'], hd_p['dec'])
	polyid_val_s      = tiling_mask.polyid(hd_s['plug_ra'], hd_s['plug_dec'])

	polyid_col_p = fits.Column(name="tiling_polyid",format='I', array=polyid_val_p )
	polyid_col_s = fits.Column(name="tiling_polyid",format='I', array=polyid_val_s )

	#number of observations in 0.5 arcsec
	# 0.5 arcsecond
	nobs_val_p = n.zeros_like(hd_p['ra'])
	nobs_val_s = n.zeros_like(hd_s['plug_ra'])
	id_ps = treeP.query_ball_tree(treeS, rmax)
	nobs_val_p = n.array([len(el) for el in id_ps])

	id_ss = treeS.query_ball_tree(treeS, rmax)
	nobs_val_s = n.array([len(el) for el in id_ss])

	nobs_col_p = fits.Column(name="nobs",format='I', array=nobs_val_p )
	nobs_col_s = fits.Column(name="nobs",format='I', array=nobs_val_s )

	# assign redshift and redshift error to photometric catalog
	z_val = n.ones_like(hd_p['ra'])*-1.
	zerr_val = n.ones_like(hd_p['ra'])*-1.
	print "merges redshift measurements with photometry"
	for ii, el in enumerate(id_ps):
		if len(el)==0:
			z_val[ii] = -9.
			zerr_val[ii] = -9.
		elif len(el)==1:
			z_val[ii] = hd_s['Z'][el[0]]
			zerr_val[ii] = hd_s['Z_ERR'][el[0]]
			print z_val[ii], zerr_val[ii]
		elif len(el)>=2:
			# gets the first guess
			z_0 = hd_s['Z'][el[0]]
			zerr_0 = hd_s['Z_ERR'][el[0]]
			# gets the second guess
			z_1 = hd_s['Z'][el[1]]
			zerr_1 = hd_s['Z_ERR'][el[1]]
			print "observed many times:", len(el), z_0, zerr_0, z_1, zerr_1
			# measurments agree
			if abs(z_0-z_1)<0.05 :
				z_val[ii] = z_0
				zerr_val[ii] = zerr_0
			else :
				id_zq = n.argmin([zerr_0, zerr_1])
				z_val[ii] = [z_0, z_1][id_zq]
				zerr_val[ii] = [zerr_0, zerr_1][id_zq]


	z_col_p = fits.Column(name="Z",format='D', array=z_val )
	zerr_col_p = fits.Column(name="Z_ERR",format='D', array=zerr_val )

	hd_p_columns = hd_p.columns + nobs_col_p + polyid_col_p + z_col_p + zerr_col_p  
	hdu = fits.BinTableHDU.from_columns(hd_p_columns)
	if os.path.isfile(output_phot_cat):
		os.remove(output_phot_cat)

	hdu.writeto(output_phot_cat)



spec_file=path_2_specObjAll#, 
geometry_file = path_2_mask_ply#)

phot_file = path_2_target_AGN_XMMSL_eboss#, 
output_phot_cat = phot_file[:-8]+"_merged.fits"
construct_summary_files_data(spec_file, phot_file, geometry_file, output_phot_cat)

phot_file = path_2_target_AGN_RASS_eboss  
output_phot_cat = phot_file[:-8]+"_merged.fits"
construct_summary_files_data(spec_file, phot_file, geometry_file, output_phot_cat)

phot_file = path_2_target_AGN_RASS_sequels 
output_phot_cat = phot_file[:-8]+"_merged.fits"
construct_summary_files_data(spec_file, phot_file, geometry_file, output_phot_cat)


def construct_random_files(spec_file, geometry_file, output_cat, FACTOR =40 ):
	print os.path.isfile(spec_file), spec_file
	print os.path.isfile(geometry_file), geometry_file
	hd_s = fits.open(spec_file)[1].data
	Nspec = len(hd_s)
	print "N targeted in chunk = ", Nspec

	mng = mangle.Mangle(geometry_file)
	raR, decR = mng.genrand(Nspec*FACTOR)
	polyid_randoms = mng.polyid(raR, decR )

	ids_with_data = n.array(list(set(hd_s['tiling_polyid'][(hd_s['tiling_polyid']>0)])))

	mask = n.in1d(polyid_randoms, ids_with_data)
	print len(raR[mask])*1./Nspec
	ra_col = fits.Column(name="RA",format='D', array=raR[mask] )
	dec_col = fits.Column(name="DEC",format='D', array=decR[mask] )
	polyid_col_r = fits.Column(name="tiling_polyid",format='I', array=polyid_randoms[mask] )

	hd_r_columns  = fits.ColDefs([ ra_col, dec_col, polyid_col_r ])  
	hdu = fits.BinTableHDU.from_columns(hd_r_columns)

	if os.path.isfile(output_cat):
		os.remove(output_cat)
	hdu.writeto(output_cat)

# generate random points in the mask

geometry_file = path_2_mask_ply

spec_file = path_2_target_AGN_XMMSL_eboss[:-8]+"_merged.fits"
output_cat = spec_file[:-5]+"_random.fits"
construct_random_files(spec_file, geometry_file, output_cat, FACTOR = 200. )

spec_file = path_2_target_AGN_RASS_eboss[:-8]+"_merged.fits"
output_cat = spec_file[:-5]+"_random.fits"
construct_random_files(spec_file, geometry_file, output_cat, FACTOR = 40. )

spec_file = path_2_target_AGN_RASS_sequels[:-8]+"_merged.fits"
output_cat = spec_file[:-5]+"_random.fits"
construct_random_files(spec_file, geometry_file, output_cat, FACTOR = 400. )


def create_TSR_catalogs_data(spec_file, completeness_per_polyid_file):
	output_phot_cat = spec_file[:-5]+"_TSR.fits"  
	hd_obs     = fits.open(spec_file)[1].data

	noFiber = (hd_obs['Z']==-9.)
	withFiber = (hd_obs['Z']>-9.)

	# compute TSR per polyID
	tsr_val_p = n.ones_like(hd_obs['ra'])*-1.
	ssr_val_p = n.ones_like(hd_obs['ra'])*-1.
	nin_val_p = n.zeros_like(hd_obs['ra'])

	all_pols = n.array(list(set(hd_obs['tiling_polyid'])))
	all_pols.sort()
	all_polys = all_pols[1:]
	tsr_poly = n.ones_like(all_polys)*-1.
	ssr_poly = n.ones_like(all_polys)*-1.
	nin_poly = n.zeros_like(all_polys)

	ns = lambda sel : float(len(sel.nonzero()[0]))
	good_z = (hd_obs['Z_ERR']>0) & (hd_obs['Z']>hd_obs['Z_ERR'])

	def get_TSR_SSR(polyID):
		inside=(hd_obs['tiling_polyid']==polyID)
		Nin = ns(inside)
		Nin_noF = ns(inside & noFiber)
		Nin_wF = ns(inside & withFiber)
		Nin_wF_gz = ns(inside & withFiber & good_z )
		return Nin, Nin_noF, Nin_wF, Nin_wF_gz, (hd_obs['tiling_polyid']==polyID)
		
	for jj, polyID in enumerate(all_polys):
		Nin, Nin_noF, Nin_wF, Nin_wF_gz, sel = get_TSR_SSR(polyID)
		#print Nin, Nin_noF, Nin_wF, Nin_wF_gz, sel
		nin_val_p[sel] = n.ones_like(tsr_val_p[sel])*Nin
		nin_poly[jj] = Nin
		if Nin>0:
			tsr_val_p[sel]=n.ones_like(tsr_val_p[sel])*Nin_wF/Nin
			tsr_poly[jj] = Nin_wF/Nin		
		if Nin>0 and Nin_wF>0:
			ssr_val_p[sel & good_z]=n.ones_like(tsr_val_p[sel & good_z])*Nin_wF_gz/Nin_wF
			ssr_poly[jj] = Nin_wF_gz/Nin_wF


	w = tsr_poly*ssr_poly

	tsr_col_p = fits.Column(name="TSR",format='D', array=tsr_val_p )
	ssr_col_p = fits.Column(name="SSR",format='D', array=ssr_val_p )
	nin_col_p = fits.Column(name="nin",format='D', array=nin_val_p )

	hd_obs_columns = hd_obs.columns + tsr_col_p + ssr_col_p + nin_col_p

	hdu = fits.BinTableHDU.from_columns(hd_obs_columns)
	if os.path.isfile(output_phot_cat):
		os.remove(output_phot_cat)

	hdu.writeto(output_phot_cat)
			
	n.savetxt(completeness_per_polyid_file, n.transpose([all_polys, nin_poly, tsr_poly, ssr_poly]), header='polyID Nobjects_inPoly TSR SSR ', fmt=['%i', '%i', '%1.5f', '%1.5f'])


spec_file = path_2_target_AGN_XMMSL_eboss[:-8]+"_merged.fits"
completeness_per_polyid_file = spec_file[:-5]+"_completeness.txt"
create_TSR_catalogs_data(spec_file, completeness_per_polyid_file)

spec_file = path_2_target_AGN_RASS_eboss[:-8]+"_merged.fits"
completeness_per_polyid_file = spec_file[:-5]+"_completeness.txt"
create_TSR_catalogs_data(spec_file, completeness_per_polyid_file)

spec_file = path_2_target_AGN_RASS_sequels[:-8]+"_merged.fits"
completeness_per_polyid_file = spec_file[:-5]+"_completeness.txt"
create_TSR_catalogs_data(spec_file, completeness_per_polyid_file)


