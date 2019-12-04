import astropy.io.fits as fits
import time

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
import numpy as n
import os
import sys

# Balmer lines

B_wave = lambda upper_level : 3645.0682 * upper_level**2 / (upper_level**2 - 2**2)

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

#em_lines = em_lines_MS
em_lines = em_lines_VdB
#abs_lines = abs_lines_MS
abs_lines = abs_lines_VdB

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
path_2_spec_dir = os.path.join(os.environ['GIT_MAKESAMPLE'], 'data/stacks')
fig_dir = path_2_spec_dir # os.path.join(os.environ['GIT_MAKESAMPLE'], 'figures/agn/figures_VAC')
#path_2_fly_dir =  os.path.join(os.environ['HOME'], 'SDSS/stacks/SPIDERS_C_GAL', 'firefly/')

name_dict = { 
	'bin1' : 'A1',
	'bin4' : 'A2', 
	'bin7' : 'A3', 
	'bin9' : 'A4', 
	'bin10': 'A5',
	'bin2' : 'B1', 
	'bin5' : 'B2',
	'bin8' : 'B3',
	'bin3' : 'B1+',
	'bin6' : 'B2+'}


path_2_stacklist_dir = os.path.join(os.environ['GIT_MAKESAMPLE'], 'data/stackLists')
#bin10.txt  bin1.txt  bin2.txt  bin3.txt  bin4.txt  bin5.txt  bin6.txt  bin7.txt  bin8.txt  bin9.txt


seq1 = [1, 2, 3]
seq2 = [1,4,7,9,10]

file_list_1 = n.array([ os.path.join( path_2_spec_dir ,'bin'+str(id_spec)+'.stack') for id_spec in seq1 ])
file_list_2 = n.array([ os.path.join( path_2_spec_dir ,'bin'+str(id_spec)+'.stack') for id_spec in seq2 ])

ascii_list_1 = n.array([ n.loadtxt(os.path.join( path_2_stacklist_dir ,'bin'+str(id_spec)+'.txt'), unpack=True)[3] for id_spec in seq1 ])
ascii_list_2 = n.array([ n.loadtxt(os.path.join( path_2_stacklist_dir ,'bin'+str(id_spec)+'.txt'), unpack=True)[3] for id_spec in seq2 ])

baseNames_1 = n.array([ name_dict[os.path.basename(bn)[:-6]] for bn in file_list_1 ])
baseNames_2 = n.array([ name_dict[os.path.basename(bn)[:-6]] for bn in file_list_2 ])


file_list = n.array( glob.glob( os.path.join( path_2_spec_dir ,'*.stack')))
file_list.sort()
print('seq1', file_list_1)
print('seq2', file_list_2)

def plot_NZ(ascii_list, baseName, path_2_figure):
	zbins = n.arange(0., 1.2, 0.1)
	fig=p.figure(1, (5.5, 5.5))
	for el, name in zip(ascii_list[::-1], baseName[::-1]):
		p.hist(el, bins = zbins, histtype='step', lw=2, label=name+ ', '+str(n.round(n.mean(el),2)))
	p.xlabel('redshift') 
	#p.ylabel(r'N')
	p.grid()
	p.legend(loc=1, fontsize=14)
	#p.tight_layout()
	p.savefig(path_2_figure)
	p.clf()

def plot_stacks(file_list_i, baseNames_i, path_2_figure='x.png', wmin=2000, wmax=9000, title_str=' '):
	print(path_2_figure, wmin, wmax)
	fig=p.figure(0, (7.5, 4.5), frameon=False )
	p.axes([0.15, 0.15, 0.82, 0.75])
	#p.title(baseNames_i[0].split('_')[0]) 
	p.xlabel('wavelength [Angstrom, rest frame]') 
	p.ylabel(r'Flux [$f_\lambda$]')
	p.grid()
	min_y_data = []
	max_y_data = []
	for jj, (path_2_spec, baseName) in enumerate(zip(file_list_i[::-1], baseNames_i[::-1])):
		d = fits.open(path_2_spec)
		sel = (d[1].data['NspectraPerPixel'] > 0.5*n.max(d[1].data['NspectraPerPixel'])) & (d[1].data['medianStack']>0)& (d[1].data['wavelength']>wmin)&(d[1].data['wavelength']<wmax)
		x_data = d[1].data['wavelength'][sel]
		y_data_i = d[1].data['medianStack'][sel]
		y_data = y_data_i/n.median(y_data_i)#/1.3**jj
		y_err = d[1].data['jackknifStackErrors'][sel]
		N_spec = d[1].data['NspectraPerPixel'][sel]
		p.plot(x_data, y_data, lw=1.5, label=baseName + ', N='+str(int(n.median(N_spec))) )
		min_y_data.append(n.min(y_data) )
		max_y_data.append(n.max(y_data) )

	min_y_data = n.array( min_y_data )
	max_y_data = n.array( max_y_data )

	em_lines.sel2 = (em_lines.sel) & (em_lines.WL>wmin) & (em_lines.WL<wmax)
	abs_lines.sel2 = (abs_lines.sel) & (abs_lines.WL>wmin) & (abs_lines.WL<wmax)

	for name, xx  in zip(em_lines.ID[em_lines.sel2], em_lines.WL[em_lines.sel2] ):
		y_id = n.searchsorted( x_data, xx)
		y_position = n.max(y_data[n.max([y_id-5,0]):n.min([y_id+5,len(x_data)])])
		print(xx, name, y_id, y_position)
		p.text(x=xx, y=y_position*1.03, s=name, withdash=True, rotation=90, fontsize=9)#, color='red')
		#p.plot(xx, y_position*1.01, marker='|', color='k')# s=name, withdash=True, rotation=90, fontsize=9)#, color='red')

	for name, xx in zip(abs_lines.ID[abs_lines.sel2], abs_lines.WL[abs_lines.sel2] ):
		y_id = n.searchsorted( x_data, xx)
		y_position = n.min(y_data[n.max([y_id-5,0]):n.min([y_id+5,len(x_data)])])
		print(xx, name, y_id, y_position)
		p.text(x=xx, y = y_position-0.1, s=name, withdash=True, rotation=90, fontsize=9)#, color='blue')
		#p.plot(xx, y_position-0.08, marker='|', color='k')# s=name, withdash=True, rotation=90, fontsize=9)#, color='red')
	
	p.axvline(6564.5377, lw=2, ls='dashed', color='k') # B_wave(3)
	p.axvline(4861.3615, lw=2, ls='dashed', color='k') # B_wave(4)
	p.axvline(4340.462 , lw=2, ls='dashed', color='k') # B_wave(5)
	p.axvline(4101.74  , lw=2, ls='dashed', color='k') # B_wave(6)
	p.axvline(3970.072 , lw=2, ls='dashed', color='k') # B_wave(7)
	p.xlim((wmin, wmax))
	p.ylim((n.min(min_y_data)-0.1, n.max(max_y_data) + 0.1))
	#p.yscale('log')
	p.title(title_str)
	p.legend(loc=2, fontsize=14)
	p.tight_layout()
	p.savefig(path_2_figure)
	p.clf()

file_list = file_list_1
prefix = 'seq1_'
baseNames = baseNames_1
ascii_list = ascii_list_1

plot_NZ(ascii_list, baseNames, os.path.join(fig_dir, prefix+'NZ.png'))

wmin = 3000
wmax = 7000
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax )

#wmin = 3900
#wmax = 4050
#plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'Hepsilon_EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax, title_str=r'H$\epsilon$ @ 3970$\AA$' )

wmin = 4000
wmax = 4200
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'Hdelta_EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax, title_str=r'H$\delta$ @ 4101$\AA$' )

wmin = 4300
wmax = 4400
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'Hgamma_EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax, title_str=r'H$\gamma$ @ 4340$\AA$' )

wmin = 4800
wmax = 5050
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'Hbeta_EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax, title_str=r'H$\beta$ @ 4861$\AA$' )

wmin = 6450
wmax = 6650
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'Halpha_EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax, title_str=r'H$\alpha$ @ 6564$\AA$' )

wrange = n.arange(3100, 8000, 1000)
for wmin, wmax in zip(wrange[:-1], wrange[1:]):
	plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax )
	
wrange = n.arange(3100, 7500, 500)
for wmin, wmax in zip(wrange[:-1], wrange[1:]):
	plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax )


file_list = file_list_2
prefix = 'seq2_'
baseNames = baseNames_2
ascii_list = ascii_list_2

plot_NZ(ascii_list, baseNames, os.path.join(fig_dir, prefix+'NZ.png'))

wmin = 3000
wmax = 7000
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax )

wmin = 4000
wmax = 4200
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'Hdelta_EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax, title_str=r'H$\epsilon$ @ 3970$\AA$' )

wmin = 4000
wmax = 4200
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'Hdelta_EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax, title_str=r'H$\delta$ @ 4101$\AA$' )

wmin = 4300
wmax = 4400
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'Hgamma_EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax, title_str=r'H$\gamma$ @ 4340$\AA$' )

wmin = 4800
wmax = 5050
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'Hbeta_EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax, title_str=r'H$\beta$ @ 4861$\AA$' )

wmin = 6450
wmax = 6650
plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'Halpha_EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax, title_str=r'H$\alpha$ @ 6564$\AA$' )

wrange = n.arange(3100, 8000, 1000)
for wmin, wmax in zip(wrange[:-1], wrange[1:]):
	plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax )
	
wrange = n.arange(3100, 7500, 500)
for wmin, wmax in zip(wrange[:-1], wrange[1:]):
	plot_stacks(file_list, baseNames, os.path.join(fig_dir, prefix+'EVstacks_wrange_'+str(wmin)+'_'+str(wmax)+'.png'), wmin, wmax )


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
		#p.xlim((x_data.min()-1, x_data.max()+1)) 
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
	wrange = n.arange(3100, 10000,1000)
	for wmin, wmax in zip(wrange[:-1], wrange[1:]):
		#wmin=3000
		#wmax=4000
		path_2_figure = os.path.join(fig_dir, os.path.basename(path_2_spec)[:-6]+'_wrange_'+str(wmin)+'_'+str(wmax)+'.png')
		plot_single_stack_with_lines(path_2_spec, wmin, wmax, path_2_figure, em_lines, abs_lines)

	