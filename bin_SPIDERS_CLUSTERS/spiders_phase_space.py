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


cat = fits.open(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'CODEX-DR14-MergedSpectroscopicCatalogue_2018-02-05.fits'))[1].data

all_clus_id = n.array(list(set(cat['CLUS_ID'])))
all_clus_id.sort()

#spm26 = fits.open(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'CODEX-DR14-MergedSpectroscopicCatalogue_2018-02-05-spm26.fits'))[1].data

#spmV5 = fits.open(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'CODEX-DR14-MergedSpectroscopicCatalogue_2018-02-05-spmv5.fits'))[1].data

spmV5 = fits.open('/home/comparat/data/spm/firefly/FireflyGalaxyEbossDR14_Chabrier_MILES.fits')[1].data

spm26 = fits.open('/home/comparat/data/spm/firefly/FireflyGalaxySdss26_Chabrier_MILES.fits')[1].data

#jjj=0
#sel = ( cat['CLUS_ID'] == all_clus_id[jjj] )
#sel26 = ( spm26['CLUS_ID'] == all_clus_id[jjj] )
#selV5 = ( spmV5['CLUS_ID'] == all_clus_id[jjj] )

#print(spm26['CLUS_ID'][sel26]
#,spmV5['CLUS_ID'][selV5]
#,cat['CLUS_ID'][sel])

#z_clus = cat['CLUZSPEC'][sel]
#dz_26 = spm26['Z'][sel26] - z_clus
#dz_V5 = spmV5['Z_NOQSO'][selV5] - z_clus


#center = coord.ICRS(ra=cat['RA_OPT'][sel]*u.degree, dec=cat['DEC_OPT'][sel]*u.degree)

#position_26 = coord.ICRS(spm26['PLUG_RA'][sel26]*u.degree, spm26['PLUG_DEC'][sel26]*u.degree)
#position_V5 = coord.ICRS(spmV5['PLUG_RA'][selV5]*u.degree, spmV5['PLUG_DEC'][selV5]*u.degree)

#sep_r200c = center.separation(position)/(cat['R200C_DEG'][sel]*u.degree)).value 



volume_rough = aa.comoving_volume(0.5)*2200.*n.pi/129600
volume = volume_rough.value

# get cluster center

# distance to center
# rescale to r200c_deg
# get the latest min(ages) of the ssp
# compute SFR

# now looks at individual galaxies
# and gets the highest SFR for each galaxy
# youngest age

spm = spm26

#for cc in cat:
cc = cat[0]

def get_props(cc, spm, z_name_spm):
	# initialize the properties to measure
	highest_sfrs = []
	youngest_ages = []
	sep_r200c = []
	delta_z = []
	stellar_masses = []
	is_member = []
	# gets the cluster properties
	z_clus = cc['CLUZSPEC']
	center = coord.ICRS(ra=cc['RA_OPT']*u.degree, dec=cc['DEC_OPT']*u.degree)
	gal = (cc['ALLFIBERID']>0)
	Ngal = len(gal.nonzero()[0])
	# loops over the galaxy members
	for id_cc, (pla, mjd, fib, is_M) in enumerate(zip(cc['ALLPLATE'][gal], cc['ALLMJD'][gal], cc['ALLFIBERID'][gal], cc['ISMEMBER'][gal])):
		# selects members
		sel = (spm['PLATE']==pla) & (spm['MJD']==mjd) & (spm['FIBERID']==fib)
		if  len(sel.nonzero()[0])>0 :
			delta_z.append( (spm[z_name_spm][sel][0] - z_clus)*3e5/cc['CLUVDISPBEST'] )
			stellar_masses.append( spm['Chabrier_MILES_stellar_mass'][sel][0] )
			stellar_ages.append( spm['Chabrier_MILES_age_massW'][sel][0] )
			is_member.append(is_M)
			position = coord.ICRS(cc['ALLRA'][id_cc]*u.degree, cc['ALLDEC'][id_cc]*u.degree)
			sep_r200c.append( (center.separation(position)/(cc['R200C_DEG']*u.degree)).value )
			#n_cp = spm['Chabrier_MILES_nComponentsSSP'][sel].astype('int')[0]
			#if n_cp > 0 :
				#all_ages = n.array([ spm['Chabrier_MILES_age_ssp_'+str(ii)][sel][0] for ii in n.arange(n_cp) ])
				#all_masses = n.array([ spm['Chabrier_MILES_stellar_mass_ssp_'+str(ii)][sel][0] for ii in n.arange(n_cp) ])
				#sfr_inst = all_masses / all_ages
				#youngest_ages.append(n.min(all_ages))
				#highest_sfrs.append(n.max(sfr_inst))
	#return n.array(delta_z), n.array(highest_sfrs), n.array(youngest_ages), n.array(sep_r200c), n.array(stellar_masses), n.array(is_member)				
	return n.array(delta_z), n.array(sep_r200c), n.array(stellar_masses), n.array(is_member)

#delta_z, highest_sfrs, youngest_ages, sep_r200c, stellar_mass , is_mem = get_props(cc, spm26, 'Z')
t0 = time.time()
DATA = [] 
for cc in cat[20:35]:
	print(cc['CLUS_ID'], time.time()-t0)
	#delta_z, highest_sfrs, youngest_ages, sep_r200c, stellar_mass , is_mem = get_props(cc, spmV5, 'Z_NOQSO')
	delta_z, sep_r200c, stellar_mass , is_mem = get_props(cc, spmV5, 'Z_NOQSO')
	mem = (is_mem==1)
	if len(delta_z[mem]>0):
		#DATA.append([delta_z[mem], highest_sfrs[mem], youngest_ages[mem], sep_r200c[mem], stellar_mass[mem]])
		DATA.append([delta_z[mem], sep_r200c[mem], stellar_mass[mem]])

p.figure(1, (5,5))
p.title('SPIDERS')
for el in DATA:
	p.scatter(el[1], el[0], s=30, c=n.log10(el[2]), edgecolors='face')

p.xlabel('r/r200c')
p.ylabel(r'c$\Delta$ z (galaxy - cluster) / $\sigma_v$')
cb=p.colorbar(shrink=0.7)
cb.set_label('stellar mass')
#p.xscale('log')
#p.yscale('log')
p.xlim((0.0,1.5))
p.ylim((-4.,5.))
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'phase-space.png'))
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

