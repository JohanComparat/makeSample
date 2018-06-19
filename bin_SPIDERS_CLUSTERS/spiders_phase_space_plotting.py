import astropy.cosmology as co
aa=co.Planck15
import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import angles
#import AngularSeparation
from astropy import coordinates as coord
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import numpy as n

import os
import sys

import ClusterScalingRelations as clsr
from scipy.interpolate import interp1d
import StellarMass as sm
smhmr = sm.StellarMass()
scl = clsr.ClusterScalingRelations_Mantz2016()

dz_1, radius_1, mass_1, age_1 = n.loadtxt(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'phase-space_26.dat'))

dz_2, radius_2, mass_2, age_2 = n.loadtxt(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'phase-space_v5.dat'))

dz_i     = n.hstack((dz_1, dz_2))
radius_i = n.hstack((radius_1, radius_2))
mass_i   = n.hstack((mass_1, mass_2))
age_i    = n.hstack((age_1, age_2))

ok = (radius_i>0 )&( dz_i > -10)&(dz_i<10)&( mass_i>0)&( age_i>0)

radius = radius_i[ok]
dz = dz_i[ok]
mass = mass_i[ok]
age = age_i[ok]

z_bins = n.arange(-3., 3., 0.2)
r_bins = n.hstack(( 1e-6, 1e-2, n.arange(0., 1.5, 0.04)[1:] ))
XX, YY = n.meshgrid((r_bins[:-1]+r_bins[1:])*0.5, (z_bins[:-1]+z_bins[1:])*0.5)

p.figure(0, (5.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
HH = n.histogram2d(radius, dz, bins=[r_bins, z_bins])[0].T
p.scatter(XX[HH>1], YY[HH>1], c=HH[HH>1], s=20, edgecolors='none', marker='s' )

p.xlabel('r/r200c')
p.ylabel(r'c$\Delta$ z (galaxy - cluster) / $\sigma_v$')
cb=p.colorbar(shrink=0.7)
cb.set_label('Number')
p.xscale('log')
#p.yscale('log')
p.xlim((0.0,1.5))
p.ylim((-3.,3.))
p.title('all members')
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'phase-space.png'))
p.clf()

med_mass = n.median(mass)

s_low = (mass < med_mass)
s_high = (mass >= med_mass)

p.figure(0, (5.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
HH_low = n.histogram2d(radius[s_low], dz[s_low], bins=[r_bins, z_bins])[0].T
HH_high = n.histogram2d(radius[s_high], dz[s_high], bins=[r_bins, z_bins])[0].T

diff = HH_high-HH_low

p.scatter(XX[HH>1], YY[HH>1], c=diff[HH>1], s=30, edgecolors='none', marker='s', cmap = 'gnuplot' )

p.xlabel('r/r200c')
p.ylabel(r'c$\Delta$ z (galaxy - cluster) / $\sigma_v$')
cb=p.colorbar(shrink=0.7)
cb.set_label('N(more massive)-N(less massive)')
p.xscale('log')
#p.yscale('log')
#p.xlim((0.0,1.5))
p.xlim((0.001,2.5))
p.ylim((-3.,3.))
p.title('all members')
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'phase-space-mass-diff.png'))
p.clf()


med_age = n.median(age)

s_low = (age < med_age)
s_high = (age >= med_age)

p.figure(0, (5.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
HH_low = n.histogram2d(radius[s_low], dz[s_low], bins=[r_bins, z_bins])[0].T
HH_high = n.histogram2d(radius[s_high], dz[s_high], bins=[r_bins, z_bins])[0].T

diff = HH_high-HH_low

p.scatter(XX[HH>1], YY[HH>1], c=diff[HH>1], s=30, edgecolors='none', marker='s', cmap='gnuplot' )

p.xlabel('r/r200c')
p.ylabel(r'c$\Delta$ z (galaxy - cluster) / $\sigma_v$')
cb=p.colorbar(shrink=0.7)
cb.set_label('N(older)-N(younger)')
p.xscale('log')
#p.yscale('log')
#p.xlim((0.0,1.5))
p.xlim((0.001,2.5))
p.ylim((-3.,3.))
p.title('all members')
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'phase-space-age-diff.png'))
p.clf()


p.figure(1, (5,5))
p.title('SPIDERS')
p.scatter(radius, dz, s=3, c=n.log10(mass), edgecolors='face')
p.xlabel('r/r200c')
p.ylabel(r'c$\Delta$ z (galaxy - cluster) / $\sigma_v$')
cb=p.colorbar(shrink=0.7)
cb.set_label('mass')
p.xscale('log')
#p.yscale('log')
p.xlim((0.000001,2.5))
p.ylim((-3.,3.))
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'phase-space-mass.png'))
p.clf()

p.figure(1, (5,5))
p.title('SPIDERS')
p.scatter(radius, dz, s=3, c=n.log10(age), edgecolors='face')
p.xlabel('r/r200c')
p.ylabel(r'c$\Delta$ z (galaxy - cluster) / $\sigma_v$')
cb=p.colorbar(shrink=0.7)
cb.set_label('age')
p.xscale('log')
#p.yscale('log')
p.xlim((0.000001,2.5))
p.ylim((-3.,3.))
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'phase-space-age.png'))
p.clf()


sys.exit()


dx = ( n.max(sep_r200c) - n.min(sep_r200c) ) /3.
r_b = n.arange(n.min(sep_r200c), n.max(sep_r200c) + dx, dx)

p.figure(1, (5,5))
for ii,bb in enumerate(r_b[:-1]):
	sub = (sep_r200c>bb)&(sep_r200c<r_b[ii+1])
	p.hist(highest_sfrs[sub], label=str(n.round(bb,3))+"<"+str(n.round(r_b[ii+1],3)), cumulative=True, normed=True, histtype='step')

p.ylabel('normed cumulative distribution')
p.xlabel('SFR [Msun/yr]')
p.xscale('log')
p.ylim((-0.01, 1.01))
p.grid()
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'disteance-2-center-SFR-histograms.png'))
p.clf()


p.figure(1, (5,5))
p.title('SPIDERS')
p.plot(sep_r200c, youngest_ages, 'r+')
p.xlabel('r/r200c')
p.ylabel('age [yr]')
p.xscale('log')
p.yscale('log')
p.xlim((0.1,5))
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'disteance-2-center-AGE.png'))
p.clf()


p.figure(1, (5,5))
p.title('SPIDERS DR14 galaxies')
p.plot(spm['Z'], spm["Chabrier_MILES_stellar_mass"], 'b,', label='targets')
p.plot(z, y, 'r,', label='cluster members')
p.xlabel('redshift')
p.ylabel('stellar mass [Msun]')
#p.xscale('log')
p.yscale('log')
p.xlim((0,0.7))
p.ylim((1e9,1e12))
p.grid()
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster',  'redshift-mass.png'))
p.clf()

logm2x = n.hstack((m2x))	
bins=n.arange(-7, 0.5, 0.1)

basis = (n.isnan(logm2x)==False)&(logm2x != -n.inf)&(logm2x != n.inf)

arbitrary_factor =5.

p.figure(1, (5,5))
ok = (basis)&(x>1e44)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])
p.plot((bins[1:]+bins[:-1])/2., n.log10(out/arbitrary_factor), label='LX>44')

ok = (basis)&(x>10**44.5)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])
p.plot((bins[1:]+bins[:-1])/2., n.log10(out/arbitrary_factor), label='LX>44.5')

ok = (basis)&(x>1e45)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])
p.plot((bins[1:]+bins[:-1])/2., n.log10(out/arbitrary_factor), label='LX>45')

ok = (basis)&(m200c>10**14)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])
p.plot((bins[1:]+bins[:-1])/2., n.log10(out/arbitrary_factor), label='M200c>14', ls='dashed')

ok = (basis)&(m200c>10**15)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])
p.plot((bins[1:]+bins[:-1])/2., n.log10(out/arbitrary_factor), label='M200c>15', ls='dashed')

xs = n.arange(-7, 0.01, 0.01)
logfsat= lambda logxi, a, b, logN0, exponent : n.log10( 10**logN0 * (10**logxi)**a)# * n.e**(-b*(10**logxi)**exponent))

p.plot(xs, logfsat(xs, -0.81, 5.81, -2.25, -2.54), label='-0.81')
p.plot(xs, logfsat(xs, -0.18, 5.81, -1.2, -.54), label='-0.18')

p.xlabel('log10(SMHMR(stellar mass) / HaloMass(Lx ray))')
p.ylabel('histogram')
#p.xscale('log')
#p.yscale('log')
p.ylim((-1.5, 0.5))
p.xlim((-4,0))
p.grid()
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'LX-mass-histogram.png'))
p.clf()

