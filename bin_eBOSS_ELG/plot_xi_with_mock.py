import glob
import sys
import os
import numpy as n
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

models = n.array(glob.glob('/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_C4/*EBOSS*.2pcf'))
models.sort()
cmass = n.loadtxt('/data36s/comparat/SDSS/LSS/LRG/cmass_boss/DR12v5_CMASS_North_rd0.monopole.2pcf', unpack=True)

topdir = "/home/comparat/data1/SDSS/ELG/clustering/" 

REF = n.loadtxt(topdir+"ebosselg.JC28Nov2017_2pcf", unpack=True)

p.figure(1, (5,5))

p.plot(REF[0], REF[0]**2*REF[1], label='eBOSS ELG')
DATA = n.loadtxt(topdir+"ebosselg.JC28Nov2017_2pcf.mask_O2O3Hb_th_1.5.txt", unpack=True)
p.plot(DATA[0], DATA[0]**2*DATA[1], label='O2O3Hb_th_1.5')
DATA = n.loadtxt(topdir+"ebosselg.JC28Nov2017_2pcf.mask_O2O3_th_1.5.txt", unpack=True)
p.plot(DATA[0], DATA[0]**2*DATA[1], label='O2O3_th_1.5')
DATA = n.loadtxt(topdir+"ebosselg.JC28Nov2017_2pcf.mask_O2_th_1.5.txt", unpack=True)
p.plot(DATA[0], DATA[0]**2*DATA[1], label='O2_th_1.5')

#for p2m in models[1::2]:
	#DATA = n.loadtxt(p2m, unpack=True)
	#p.plot(DATA[0], DATA[0]**2*DATA[1], lw=0.8, label=os.path.basename(p2m)[:-5] )
p2m=models[3]
DATA = n.loadtxt(p2m, unpack=True)
p.plot(DATA[0], DATA[0]**2*DATA[1], lw=0.8, label=os.path.basename(p2m)[:-5] )
p2m=models[5]
DATA = n.loadtxt(p2m, unpack=True)
p.plot(DATA[0], DATA[0]**2*DATA[1], lw=0.8, label=os.path.basename(p2m)[:-5] )

#p.xscale('log')
#p.yscale('log')
p.xlim((0.2, 160))
p.ylim((-20,60))
p.xlabel('s [Mpc/h]')
p.ylabel(r'$s^2.\xi(s)$')
p.grid()
p.legend(frameon=False, loc=0, fontsize=10)
p.savefig('/home/comparat/wwwDir/eRoMok/elg_clustering.png')
p.clf()

from scipy.interpolate import interp1d

p.figure(1, (5,5))

#DATA = n.loadtxt(topdir+"ebosselg.JC28Nov2017_2pcf.mask_O2O3Hb_th_1.5.txt", unpack=True)
dd = interp1d(REF[0], REF[1])
p.plot(REF[0], 1+REF[2]**(-0.5), 'k--', lw=0.5, label='error' )
p.plot(REF[0], 1-REF[2]**(-0.5), 'k--', lw=0.5 )

for p2m in models[1::2]:
  DATA = n.loadtxt(p2m, unpack=True)
  sel = (DATA[0]>0.2)&(DATA[0]<80)
  y_data = dd(DATA[0][sel])
  p.plot(DATA[0][sel], y_data/DATA[1][sel], lw=0.8, label=os.path.basename(p2m)[:-5] )


#p.xscale('log')
#p.yscale('log')
p.xlim((0.2, 60))
p.ylim((0.5, 1.5))
p.xlabel('s [Mpc/h]')
p.ylabel(r'$\xi(s)$ data / model')
p.grid()
p.legend(frameon=False, loc=0, fontsize=10)
p.savefig('/home/comparat/wwwDir/eRoMok/elg_clustering_ratio.png')
p.clf()


