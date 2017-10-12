import astropy.io.fits as fits
import os
import sys
from os.path import join
import numpy as n
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

d5 = n.loadtxt("/data36s/comparat/SDSS/targets/ELG/DR5.2pcf"  , unpack=True)
d52 = n.loadtxt("/data36s/comparat/SDSS/targets/ELG/DR5_2.2pcf", unpack=True)
f16 = n.loadtxt("/data36s/comparat/SDSS/targets/ELG/favole_angular.2pcf"  , unpack=True)


p.figure(1, (5,5))
p.axes([0.17,0.17,0.78,0.78])
p.plot(d5[0],d5[1],label='dr5-eboss')
p.plot(d52[0],d52[1],label='dr5-eboss2')
p.plot(f16[0],f16[1],label='mock elg')
p.xlabel(r'$\theta$ [deg]')
p.ylabel(r'$w(\theta)$')
p.xscale('log')
p.yscale('log')
p.xlim((0.001,0.5))
p.ylim((0.01, 2))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig(join(os.environ['HOME'], 'wwwDir/sdss/elg/wtheta_targets.png'))
p.clf()
