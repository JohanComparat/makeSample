import os, sys, glob
from os.path import join
import numpy as np
import matplotlib.pyplot as p

out_dir = '/data36s/comparat/CODEX_clustering/angular_clustering/'

xcorr_list = np.array(glob.glob(os.path.join(out_dir , 'cat_spiders_masked_X_GAIA_table_*.fits.data')))
xcorr_list.sort()

valid_ids = np.array([0, 1, 2,9,10,11,12,13])

def smooth(x,window_len=11,window='hanning'):
	"""
	window can be : 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
	"""
	s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
	if window == 'flat': 
		w=np.ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')
	y=np.convolve(w/w.sum(),s,mode='valid')
	return y[(window_len/2-1):-(window_len/2)]


window_len=4

rho = 0.0287

p.figure(1, (8,5))

for el in xcorr_list[:7]:
	print(el)
	x,y=np.loadtxt(el, unpack=True)
	y1=smooth(y, window_len=window_len, window='blackman')
	lab = os.path.basename(el)[:-10].split('_')
	p.plot(x*3600.,y/rho, label=lab[-3]+'<g<'+lab[-1], lw=2)
	#p.plot(x*3600.,y1/rho, label='sm '+lab[-3]+'<g<'+lab[-1], lw=2)

p.legend(frameon=False, loc=0)
p.xlabel('separation to star [arcseconds]')
p.ylabel('relative density of targets')
p.xscale('log')
p.xlim((0.1,40))
p.ylim((0.,2.))
#p.yscale('log')
p.grid()
p.savefig(os.path.join(out_dir ,'CODEX_GAIA_faint.png'))
p.clf()

p.figure(1, (8,5))

for el in xcorr_list[7:]:
	print(el)
	x,y=np.loadtxt(el, unpack=True)
	y1=smooth(y, window_len=window_len, window='blackman')
	lab = os.path.basename(el)[:-10].split('_')
	p.plot(x*3600.,y/rho, label=lab[-3]+'<g<'+lab[-1], lw=2)
	#p.plot(x*3600.,y1/rho, label='sm '+lab[-3]+'<g<'+lab[-1], lw=2)

p.legend(frameon=False, loc=0)
p.xlabel('separation to star [arcseconds]')
p.ylabel('relative density of targets')
p.xscale('log')
p.xlim((0.1,200))
p.ylim((0.,2.))
#p.yscale('log')
p.grid()
p.savefig(os.path.join(out_dir ,'CODEX_GAIA_bright.png'))
p.clf()


os.system("cp "+out_dir+"CODEX_GAIA*.png /home/comparat/wwwDir/stuff/")


