import astropy.io.fits as fits
import time

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
import numpy as n
import os
import sys

class Lines:
    pass

path_2_abs_line = os.path.join(os.environ['GIT_MAKESAMPLE'],'data/line-list-absorption-vandenberk2001.txt')  
path_2_em_line = os.path.join(os.environ['GIT_MAKESAMPLE'],'data/line-list-emission-vandenberk2001.txt')

abs_lines_VdB = Lines()
wl, wl_er, W, W_er, Width, ID, rest_wl = n.loadtxt(path_2_abs_line, unpack=True, dtype={'names': ('wl', 'wl_er', 'W', 'W_er', 'Width', 'ID', 'rest_wl'), 'formats': ('f4', 'f4', 'f4', 'f4', 'f4', 'S10', 'f4')})
abs_lines_VdB.ID = ID.astype('str')
abs_lines_VdB.WL = rest_wl
abs_lines_VdB.sel = (abs_lines_VdB.WL>0)

em_lines_VdB = Lines()
obswl,obswlerr,wllo,wlhi,relFl,relFlerr,  W,   W_err,  Width, Skew, ID, labwl = n.loadtxt(path_2_em_line, unpack=True, dtype={'names': ('obswl', 'obswlerr', 'wllo', 'wlhi', 'relFl', 'relFlerr', '  W', '   W_err', '  Width', ' Skew', ' ID', ' labwl'), 'formats': ('f4', 'f4','f4','f4',  'f4','f4','f4',  'f4','f4','f4', 'S10', 'S10')})
em_lines_VdB.ID = ID.astype('str')
em_lines_VdB.WL = obswl
em_lines_VdB.sel = (em_lines_VdB.WL>0)


path_2_abs_line = os.path.join(os.environ['GIT_MAKESAMPLE'],'data/line-list-absorption-salvato.dat')  
path_2_em_line = os.path.join(os.environ['GIT_MAKESAMPLE'],'data/line-list-emission-salvato.dat')

abs_lines_MS = Lines()
wl, ID = n.loadtxt(path_2_abs_line, unpack=True, dtype={'names': ('wl', 'ID'), 'formats': ('f4', 'S100')})
abs_lines_MS.ID = ID.astype('str')
abs_lines_MS.WL = wl
abs_lines_MS.sel = (abs_lines_MS.WL>0)

em_lines_MS = Lines()
wl, ID = n.loadtxt(path_2_em_line, unpack=True, dtype={'names': ('wl', 'ID'), 'formats': ('f4', 'S100')})
em_lines_MS.ID = ID.astype('str')
em_lines_MS.WL = wl
em_lines_MS.sel = (em_lines_MS.WL>0)


from cycler import cycler
# 1. Setting prop cycle on default rc parameter
p.rc('lines', linewidth=1.3)
p.rc('axes', prop_cycle=cycler('color', ["#638bd9", "#e586b6", "#dd7c25", "#598664", "#9631a9", "#cff159", "#534f55", "#ab1519", "#89dbde"]) )


from scipy.stats import scoreatpercentile
from scipy.stats import norm
from scipy.interpolate import interp1d
m_bins = n.arange(-4,4,0.1)

import glob
path_2_spec_dir = os.path.join(os.environ['HOME'], 'SDSS/stacks/X_AGN')
fig_dir = os.path.join(os.environ['GIT_MAKESAMPLE'], 'figures/agn/figures_VAC')
#path_2_fly_dir =  os.path.join(os.environ['HOME'], 'SDSS/stacks/SPIDERS_C_GAL', 'firefly/')
file_list = n.array( glob.glob( os.path.join( path_2_spec_dir ,'*.stack')))
file_list.sort()
baseNames = n.array([ os.path.basename(bn)[5:-6] for bn in file_list ])
baseNamesSplit = n.array([ bn.split('_') for bn in baseNames ])
uniq_type = n.unique(baseNamesSplit.T[0])

def plot_stacks(file_list_i, baseNames_i, path_2_figure='x.png'):
	print(path_2_figure)
	p.figure(0, (12.2, 8.2), frameon=False )
	#p.axes([0.1, 0.1, 0.85, 0.78])
	#p.title(baseNames_i[0].split('_')[0]) 
	p.xlabel('wavelength [Angstrom, rest frame]') 
	p.ylabel(r'Flux [$f_\lambda$]')
	p.xlim=((800, 7500)) 
	p.grid()
	
	for jj, (path_2_spec, baseName) in enumerate(zip(file_list_i[1::3], baseNames_i[1::3])):
		d = fits.open(path_2_spec)
		sel = (d[1].data['NspectraPerPixel'] > 0.5*n.max(d[1].data['NspectraPerPixel'])) & (d[1].data['medianStack']>0)
		x_data = d[1].data['wavelength'][sel]
		y_data = d[1].data['medianStack'][sel]
		y_err = d[1].data['jackknifStackErrors'][sel]
		N_spec = d[1].data['NspectraPerPixel'][sel]
		labels = baseName.split('_')
		p.plot(x_data, y_data/2**jj, lw=1, label=str(float(labels[2])/10.)+'<z<'+str(float(labels[4])/10.) + ', N='+str(int(n.median(N_spec))) )

	#for name, xx  in zip(em_lines.ID[em_lines.sel], em_lines.WL[em_lines.sel] ):
		#print(xx, 20, name)
		#p.text(x=xx, y=20, s=name, withdash=True, rotation=90, fontsize=9)#, color='red')
		#p.plot(xx, 20, marker='|', color='k')# s=name, withdash=True, rotation=90, fontsize=9)#, color='red')

	#for name, xx in zip(abs_lines.ID[abs_lines.sel], abs_lines.WL[abs_lines.sel] ):
		#print(xx, 1, name)
		#p.text(x=xx, y=1, s=name, withdash=True, rotation=90, fontsize=9)#, color='blue')
		#p.plot(xx, 1, marker='|', color='k')# s=name, withdash=True, rotation=90, fontsize=9)#, color='red')
		
	p.yscale('log')
	p.legend(loc=0, fontsize=8, frameon=False)
	p.tight_layout()
	p.savefig(path_2_figure)
	p.clf()


for UT in uniq_type:
	sel = ( baseNamesSplit.T[0] == UT )
	plot_stacks(file_list[sel], baseNames[sel], path_2_figure=os.path.join(fig_dir, UT+'.png') )
	
sys.exit()

def plot_single_stack_with_lines(path_2_spec, wmin, wmax, path_2_figure, em_lines, abs_lines):
	print(path_2_spec, path_2_figure)
	d = fits.open(path_2_spec)
	sel = (d[1].data['NspectraPerPixel'] > 0.5*n.max(d[1].data['NspectraPerPixel'])) & (d[1].data['medianStack']>0) & (d[1].data['wavelength']>wmin)&(d[1].data['wavelength']<wmax)
	x_data = d[1].data['wavelength'][sel]
	y_data = d[1].data['medianStack'][sel]
	#print(x_data, y_data)
	y_err = d[1].data['jackknifStackErrors'][sel]
	N_spec = d[1].data['NspectraPerPixel'][sel]
	#labels = baseName.split('_')
	if len(x_data)>20:
		em_lines.sel2 = (em_lines.sel) & (em_lines.WL>wmin) & (em_lines.WL<wmax)
		abs_lines.sel2 = (abs_lines.sel) & (abs_lines.WL>wmin) & (abs_lines.WL<wmax)

		p.figure(0, (8., 3.5), frameon=False )
		# panel top left: spectrum + models
		p.axes([0.1, 0.2, 0.85, 0.78])
		p.xlabel('wavelength [Angstrom, rest frame]')
		p.ylabel(r'Flux [$f_\lambda$]')#, # 10^{-17}$ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.xlim((x_data.min()-1, x_data.max()+1)) 
		p.ylim((y_data.min()-3, y_data.max()+6)) 
		p.grid()

		p.plot(x_data, y_data, lw=2, ls='solid', color='k' )#, label=str(float(labels[2])/10.)+'<z<'+str(float(labels[4])/10.) + ', N='+str(int(n.median(N_spec))) )

		for name, xx  in zip(em_lines.ID[em_lines.sel2], em_lines.WL[em_lines.sel2] ):
			y_id = n.searchsorted( x_data, xx)
			y_position = n.max(y_data[n.max([y_id-10,0]):n.min([y_id+10,len(x_data)])])
			print(xx, name, y_id, y_position)
			p.text(x=xx, y=y_position+5, s=name, withdash=True, rotation=90, fontsize=9)#, color='red')
			p.plot(xx, y_position+1.5, marker='|', color='k')# s=name, withdash=True, rotation=90, fontsize=9)#, color='red')

		for name, xx in zip(abs_lines.ID[abs_lines.sel2], abs_lines.WL[abs_lines.sel2] ):
			y_id = n.searchsorted( x_data, xx)
			y_position = n.min(y_data[n.max([y_id-10,0]):n.min([y_id+10,len(x_data)])])
			print(xx, name, y_id, y_position)
			p.text(x=xx, y=y_position-2, s=name, withdash=True, rotation=90, fontsize=9)#, color='blue')
			p.plot(xx, y_position-1, marker='|', color='k')# s=name, withdash=True, rotation=90, fontsize=9)#, color='red')
			
		#p.yscale('log')
		p.legend(loc=0, fontsize=8, frameon=False)
		#p.tight_layout()
		p.savefig(path_2_figure)
		p.clf()

for path_2_spec in file_list:
	#em_lines = em_lines_MS
	em_lines = em_lines_VdB
	#abs_lines = abs_lines_MS
	abs_lines = abs_lines_VdB
	wrange = n.arange(100, 10000,1000)
	for wmin, wmax in zip(wrange[:-1], wrange[1:]):
		#wmin=3000
		#wmax=4000
		path_2_figure = os.path.join(fig_dir, os.path.basename(path_2_spec)[:-6]+'_wrange_'+str(wmin)+'_'+str(wmax)+'.png')
		plot_single_stack_with_lines(path_2_spec, wmin, wmax, path_2_figure, em_lines, abs_lines)

	