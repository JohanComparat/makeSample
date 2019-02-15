import os, sys, glob
from os.path import join
import numpy as np
import matplotlib.pyplot as p

out_dir = '/data36s/comparat/AGN_clustering/angular_clustering/'

xcorr_list = np.array(glob.glob(os.path.join(out_dir , '2RXS_AllWISE_catalog_paper_2017May26_X_GAIA_table_*.fits.data')))
xcorr_list.sort()


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

p.figure(1, (8,5))

for el in xcorr_list:
	print(el)
	x,y=np.loadtxt(el, unpack=True)
	y1=smooth(y, window_len=window_len, window='blackman')
	lab = os.path.basename(el)[:-10].split('_')
	p.plot(x*3600.,y, label=lab[-3]+'<g<'+lab[-1], lw=2)
	#p.plot(x*3600.,y1, label='sm '+lab[-3]+'<g<'+lab[-1], lw=2)

p.legend(frameon=False, loc=0)
p.xlabel('separation to star [arcseconds]')
p.ylabel('relative density of targets')
p.xscale('log')
p.xlim((0.1,20))
p.yscale('log')
p.grid()
p.savefig(os.path.join(out_dir ,'2RXS_GAIA.png'))
p.clf()

os.system("cp /data36s/comparat/AGN_clustering/angular_clustering/2RXS_GAIA.png /home/comparat/wwwDir/stuff/")



from numpy import *
from pylab import *

def smooth_demo():

    t=linspace(-4,4,100)
    x=sin(t)
    xn=x+randn(len(t))*0.1
    y=smooth(x)

    ws=31

    subplot(211)
    plot(ones(ws))

    windows=['flat', 'hanning', 'hamming', 'bartlett', 'blackman']

    hold(True)
    for w in windows[1:]:
        eval('plot('+w+'(ws) )')

    axis([0,30,0,1.1])

    legend(windows)
    title("The smoothing windows")
    subplot(212)
    plot(x)
    plot(xn)
    for w in windows:
        plot(smooth(xn,10,w))
    l=['original signal', 'signal with noise']
    l.extend(windows)

    legend(l)
    title("Smoothing a noisy signal")
    show()

