import numpy as n
import glob
import os, sys
import astropy.io.fits as fits
import time
t0=time.time()

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

agn_clustering_dir = '/data36s/comparat/AGN_clustering'
agn_clustering_dir = '/home/comparat/data/AGN_clustering'

catalog_dir  = os.path.join(agn_clustering_dir, 'catalogs'  )
figure_dir = os.path.join(catalog_dir, 'figures' )

if os.path.isdir(figure_dir)==False:
	os.system('mkdir -p '+figure_dir)

spec_list = n.array( glob.glob( os.path.join( os.environ['HOME'], 'SDSS/stacks/X_AGN/*zmin_00_zmax_50.asc') ) )
spec_list.sort()
list_2_plot = n.array( glob.glob( os.path.join( os.environ['HOME'], 'SDSS/stacks/X_AGN/*zmin_00_zmax_50.stack.specMatrix.dat') ) )
list_2_plot.sort()
list_stacks = n.array( glob.glob( os.path.join( os.environ['HOME'], 'SDSS/stacks/X_AGN/*zmin_00_zmax_50.stack') ) )
list_stacks.sort()

ii = 2 

path_2_stack = list_stacks[ii]
path_2_specs = list_2_plot[ii]
path_2_list  = spec_list[ii]

stack = fits.open(path_2_stack)
data = n.loadtxt(path_2_specs)
plates, mjds, fiberids, redshifts = n.loadtxt(path_2_list, unpack=True)

ids = n.argsort(redshifts)[::-1]
ll = stack[1].data['wavelength']
w_selection = (ll>2200)&(ll<8000)
w_sel = n.arange(len(ll))[(ll>2200)&(ll<8000)]
w_min = n.min(w_sel)
w_max = n.max(w_sel)

sorted_matrix = data[ids][:,w_min:w_max]
sorted_zs = redshifts[ids]
#def resample_it(IDX_j):
	##IDX_j = 122880
	#IDX_min=IDX_j
	#IDX_max=IDX_j+4096
	#IDX_str=str(IDX_min).zfill(6)+'-'+str(IDX_max).zfill(6)
	#FLUXES     = n.loadtxt(out_file+'.'+IDX_str+'.dat')
	#med_val = n.median(FLUXES, axis=1)
	#F_med = FLUXES/med_val
	#out = n.median(F_med.reshape(128,32,4096), axis=1)
	#wavelength = n.loadtxt(out_file+'.wavelength.'+IDX_str+'.dat')
	#data       = n.loadtxt(out_file+'.shapes.'+IDX_str+'.dat')
	#plates, mjds, fiberids, redshifts = n.loadtxt(out_file+'.list.'+IDX_str+'.dat')
	#out_z = n.median(redshifts.reshape(128, 32), axis=1)
	#return out, out_z, wavelength


path_2_figure = os.path.join( figure_dir, os.path.basename(path_2_specs) + '.png')
path_2_figure = os.path.join( os.path.join(os.environ['HOME'], 'wwwDir', 'stuff', os.path.basename(path_2_specs) + '.png'))

#out_0, out_z_0, wavelength = resample_it(IDX_j)
#out_1, out_z_1, wavelength = resample_it(IDX_j+4096)
#out_2, out_z_2, wavelength = resample_it(IDX_j+4096*2)
#out_3, out_z_3, wavelength = resample_it(IDX_j+4096*3)

#out = n.vstack(( out_0, out_1, out_2, out_3 ))
#out_z = n.hstack(( out_z_0, out_z_1, out_z_2, out_z_3 ))

#xx = n.searchsorted(wavelength, n.array([2800., 3728., 5007., ])*(1+redshifts[0]))
#yy = n.zeros(len(xx))
#tt = n.array(['MgII 2800', '[OII] 3728', '[OIII] 5007' ])
#out[out>100.]=100.
#out[out<0.01]=0.01
fig = p.figure(2, (14.5, 6.5) )
#p.imshow(n.log10(FLUXES/med_val), cmap='tab20', rasterized = True)
p.imshow(n.log10(sorted_matrix), cmap='gist_gray', rasterized = True, vmin=-1, vmax=2.)
p.colorbar(shrink=0.7)
idx = n.arange(0, len(ll[w_selection]), int(len(ll[w_selection])/10))
p.xticks(idx, ll[w_selection][idx].astype('int'))
idy = n.arange(0, len(sorted_zs), int(len(sorted_zs)/10))
p.yticks(idy, n.round(sorted_zs[idy],2))
p.xlabel('wavelength')
p.ylabel('redshift')
#for x,y,t in zip(xx,yy,tt):
	#p.plot(x,y,'ro', markersize=5)
	#p.text(x,y-10,t,color='b', fontsize=14, rotation=0)
p.savefig(path_2_figure)
p.clf()

