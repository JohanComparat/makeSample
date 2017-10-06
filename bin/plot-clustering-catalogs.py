#!/usr/bin/env python

import astropy.io.fits as fits
import os
import sys
from os.path import join
import numpy as n
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p

vac_dir = join( os.environ['DECALS_DATA_DIR'], "catalogs", "VAC")

def plot_ra_dec_catalog(file,aaa=1.):
	base = os.path.basename(file)[:-5]
	hd = fits.open(file)[1].data
	p.figure(0,(6,6))
	p.plot(hd['ra'], hd['dec'], 'k,', rasterized=True,alpha=aaa)
	p.xlabel(r'dec')
	p.ylabel(r'ra')
	p.grid()
	p.title(base)
	p.savefig(join(vac_dir, "plots", "plot_"+base+"_ra_dec.png"))
	p.clf()

def compare_histograms(file_data, file_random, qty):
	base = os.path.basename(file_data)[:-5]
	hd = fits.open(file_data)[1].data
	hr = fits.open(file_random)[1].data
	p.hist(hd[qty], normed=True, cumulative=True,histtype='step', bins=100, label='data')
	p.hist(hr[qty], normed=True,cumulative=True,histtype='step', bins=100, label='random')
	p.xlabel(qty)
	p.ylabel('Cumulative weighted counts')
	p.grid()
	p.title(base)
	p.legend(frameon=False, loc=0)
	p.savefig(join(vac_dir, "plots", "hist_"+qty+"_"+base+"_hist_comparison.png"))
	p.clf()

def compare_histograms_w(file_data, file_random, qty):
	base = os.path.basename(file_data)[:-5]
	hd = fits.open(file_data)[1].data
	hr = fits.open(file_random)[1].data
	p.hist(hd[qty], weights=1/hd['w'], normed=True, cumulative=True,histtype='step', bins=100, label='data')
	p.xlabel(qty)
	p.ylabel('Weighted cumulative normed counts')
	p.hist(hr[qty], normed=True,cumulative=True,histtype='step', bins=100, label='random')
	p.legend(frameon=False, loc=0)
	p.grid()
	p.title(base)
	p.savefig(join(vac_dir, "plots", "hist_"+qty+"_"+base+"_hist_comparison.png"))
	p.clf()

plot_ra_dec_catalog(join(vac_dir, "elg_240_sgc.v2.clustering.chunk21.fits"))
#plot_ra_dec_catalog(join(vac_dir, "elg_240_sgc.v2.clustering.wtheta.chunk21.fits"),aaa=0.1)
plot_ra_dec_catalog(join(vac_dir, "random-sweep.clustering.chunk21.fits"))
#plot_ra_dec_catalog(join(vac_dir, "random-sweep.clustering.wtheta.chunk21.fits"),aaa=0.1)
"""
file_data = join(vac_dir, "elg_240_sgc.v2.clustering.wtheta.chunk21.fits")
file_random = join(vac_dir, "random-sweep.clustering.wtheta.chunk21.fits")
qty = 'ra'
compare_histograms(file_data, file_random, qty)
qty = 'dec'
compare_histograms(file_data, file_random, qty)
"""
file_data = join(vac_dir, "elg_240_sgc.v2.clustering.chunk21.fits")
file_random = join(vac_dir, "random-sweep.clustering.chunk21.fits")
qty = 'ra'
compare_histograms_w(file_data, file_random, qty)
qty = 'dec'
compare_histograms_w(file_data, file_random, qty)
qty = 'z'
compare_histograms_w(file_data, file_random, qty)

plot_ra_dec_catalog(join(vac_dir, "elg_240_sgc.v2.clustering.chunk22.fits"))
#plot_ra_dec_catalog(join(vac_dir, "elg_240_sgc.v2.clustering.wtheta.chunk22.fits"))
plot_ra_dec_catalog(join(vac_dir, "random-sweep.clustering.chunk22.fits"))
#plot_ra_dec_catalog(join(vac_dir, "random-sweep.clustering.wtheta.chunk22.fits"))

"""
file_data = join(vac_dir, "elg_240_sgc.v2.clustering.wtheta.chunk22.fits")
file_random = join(vac_dir, "random-sweep.clustering.wtheta.chunk22.fits")
qty = 'ra'
compare_histograms(file_data, file_random, qty)
qty = 'dec'
compare_histograms(file_data, file_random, qty)
"""
file_data = join(vac_dir, "elg_240_sgc.v2.clustering.chunk22.fits")
file_random = join(vac_dir, "random-sweep.clustering.chunk22.fits")
qty = 'ra'
compare_histograms_w(file_data, file_random, qty)
qty = 'dec'
compare_histograms_w(file_data, file_random, qty)
qty = 'z'
compare_histograms_w(file_data, file_random, qty)
