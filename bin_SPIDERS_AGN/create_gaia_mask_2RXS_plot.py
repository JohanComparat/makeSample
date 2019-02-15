import os, sys, glob
from os.path import join
import numpy as np
import matplotlib.pyplot as p

out_dir = '/data36s/comparat/AGN_clustering/angular_clustering/'

xcorr_list = np.array(glob.glob(os.path.join(out_dir , '2RXS_AllWISE_catalog_paper_2017May26_X_GAIA_table_*.fits.data')))
xcorr_list.sort()

p.figure(1, (8,5))

for el in xcorr_list:
	print(el)
	x,y=np.loadtxt(el, unpack=True)
	lab = os.path.basename(el)[:-10].split('_')
	p.plot(x*3600.,y, label=lab[-3]+'<g<'+lab[-1], lw=2)

p.legend(frameon=False, loc=0)
p.xlabel('separation to star [arcseconds]')
p.ylabel('relative density of targets')
p.xscale('log')
p.xlim((0.5,40))
p.yscale('log')
p.grid()
p.savefig(os.path.join(out_dir ,'2RXS_GAIA.png'))
p.clf()

#cp /data36s/comparat/AGN_clustering/angular_clustering/2RXS_GAIA.png /home/comparat/wwwDir/stuff/

