import astropy.cosmology as co
aa=co.Planck15
import astropy.io.fits as fits
import astropy.units as u

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p


from astropy.coordinates import angles
#import AngularSeparation
from astropy import coordinates as coord

import numpy as n
import os

import sys

import ClusterScalingRelations as clsr
from scipy.interpolate import interp1d
import StellarMass as sm
smhmr = sm.StellarMass()
scl = clsr.ClusterScalingRelations_Mantz2016()

cat = fits.open(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'validatedclusters_catalogue_2016-07-04-DR14_version_round1-v4_Xmass-v1.fits.gz'))[1].data
spm = fits.open(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'validatedclusters_catalogue_2016-07-04-DR14_version_round1-v4_Xmass-v1_spm.fits'))[1].data

dr14 = fits.open(os.path.join(os.environ['DATA_DIR'], 'SDSS', 'specObj-dr14.fits'))

eboss = fits.open("/home/comparat/data/spm/firefly/FireflyGalaxyEbossDR14_Chabrier_MILES.fits")[1].data
eboss_ok = (n.log10(eboss['Chabrier_MILES_stellar_mass_up']) - n.log10(eboss['Chabrier_MILES_stellar_mass_low']) < 0.4 )


sdss = fits.open("/home/comparat/data/spm/firefly/FireflyGalaxySdss26_Chabrier_MILES.fits")[1].data
sdss_ok = (n.log10(sdss['Chabrier_MILES_stellar_mass_up']) - n.log10(sdss['Chabrier_MILES_stellar_mass_low']) < 0.4 )

z_2_d = interp1d(n.arange(0,1.2,0.0001), aa.comoving_distance(n.arange(0,1.2,0.0001)))
d_2_z = interp1d(aa.comoving_distance(n.arange(0,1.2,0.0001)), n.arange(0,1.2,0.0001))

volume_rough = aa.comoving_volume(0.5)*2200.*n.pi/129600
volume = volume_rough.value

N_R200 = 10.

#previous catalogs
m2x = []
lx = []
mass = []
zs = []
m200c = []
crd = []

p.figure(1, (5,5))
p.title('DR14 galaxies around clusters at 0.05<z<0.1')

for cc in cat[(cat['SCREEN_CLUZSPEC']>0.05)&(cat['SCREEN_CLUZSPEC']<0.1)]:
	#index = 0
	#cc = cat[index]
	#cc['R200C_DEG']
	center = coord.ICRS(ra=cc['RA_OPT']*u.degree, dec=cc['DEC_OPT']*u.degree)
	#gal = (spm['CLUS_ID']==cc['CLUS_ID'])
	##&(spm['Chabrier_MILES_stellar_mass']>0.)
	#member = (cc['SCREEN_ISMEMBER_W'][:len(gal.nonzero()[0])]==1)
	##&(abs(spm['Z']-cc['SCREEN_CLUZSPEC'])<0.2)
	#stellar_mass = spm['Chabrier_MILES_stellar_mass'][gal][member]
	#sm =stellar_mass[stellar_mass>0]

	# select in dr14 objects close enough
	deg_per_mpc = aa.arcsec_per_kpc_comoving(cc['SCREEN_CLUZSPEC']).value/3.6
	r200_mpc = cc['R200C_DEG'] / deg_per_mpc
	z_max = d_2_z(aa.comoving_distance(cc['SCREEN_CLUZSPEC']).value+N_R200*r200_mpc)
	z_min = d_2_z(aa.comoving_distance(cc['SCREEN_CLUZSPEC']).value-N_R200*r200_mpc)

	#around = (abs(dr14[1].data['PLUG_RA']-cc['RA_OPT'])<N_R200*cc['R200C_DEG'])&(abs(dr14[1].data['PLUG_DEC']-cc['DEC_OPT'])<N_R200*cc['R200C_DEG'])&(dr14[1].data['Z']<z_max)&(dr14[1].data['Z']>z_min)&(dr14[1].data['ZWARNING']==0)
	#len(dr14[1].data[around])

	around_sdss = (sdss_ok)&(abs(sdss['PLUG_RA']-cc['RA_OPT'])<N_R200*cc['R200C_DEG'])&(abs(sdss['PLUG_DEC']-cc['DEC_OPT'])<N_R200*cc['R200C_DEG'])&(sdss['Z']<z_max)&(sdss['Z']>z_min)&(sdss['ZWARNING']==0)
	
	print("in sdss",len(sdss[around_sdss]))
	if len(sdss[around_sdss])>0:
		ra_sdss = sdss['PLUG_RA'][around_sdss]
		dec_sdss = sdss['PLUG_DEC'][around_sdss]
		z_sdss = sdss['Z'][around_sdss]
		n_cp = sdss['Chabrier_MILES_nComponentsSSP'][around_sdss].astype('int')
		youngest_ages=[]
		highest_sfrs=[]
		for jj in n.arange(len(sdss[around_sdss])):
			if n_cp[jj] > 0 :
				all_ages = n.array([ sdss['Chabrier_MILES_age_ssp_'+str(ii)][around_sdss][jj] for ii in n.arange(n_cp[jj]) ])
				all_masses = n.array([ sdss['Chabrier_MILES_stellar_mass_ssp_'+str(ii)][around_sdss][jj] for ii in n.arange(n_cp[jj]) ])
				sfr_inst = all_masses / all_ages
				youngest_ages.append(n.min(all_ages))
				highest_sfrs.append(n.max(sfr_inst))
			
		highest_sfrs = n.array(highest_sfrs)
		youngest_ages = n.array(youngest_ages)
		age_sdss = youngest_ages # highest_sfrs
		#age_sdss = sdss['Chabrier_MILES_stellar_mass'][around_sdss]/sdss['Chabrier_MILES_age_lightW'][around_sdss]
		#age_sdss = sdss['Chabrier_MILES_age_lightW'][around_sdss]

		position_sdss = coord.ICRS(ra_sdss*u.degree, dec_sdss*u.degree)
		sep_sky_r200c_mpc = (center.separation(position_sdss)/(cc['R200C_DEG']*u.degree)).value * r200_mpc
		sep_los_r200c_mpc = abs(aa.comoving_distance(z_sdss)-aa.comoving_distance(cc['SCREEN_CLUZSPEC'])).value
		sdss_sep_3d_per_r200 = (sep_sky_r200c_mpc**2.+sep_los_r200c_mpc**2.)**(0.5)/r200_mpc

		p.plot(sdss_sep_3d_per_r200, age_sdss, 'r,', rasterized=True)

	around_eboss = (eboss_ok)&(abs(eboss['PLUG_RA']-cc['RA_OPT'])<N_R200*cc['R200C_DEG'])&(abs(eboss['PLUG_DEC']-cc['DEC_OPT'])<N_R200*cc['R200C_DEG'])&(eboss['Z']<z_max)&(eboss['Z']>z_min)&(eboss['ZWARNING']==0)
	print("in eboss", len(eboss[around_eboss]))
	if  len(eboss[around_eboss])>0:
		ra_eboss = eboss['PLUG_RA'][around_eboss]
		dec_eboss = eboss['PLUG_DEC'][around_eboss]
		z_eboss = eboss['Z'][around_eboss]
		n_cp = eboss['Chabrier_MILES_nComponentsSSP'][around_eboss].astype('int')
		youngest_ages=[]
		highest_sfrs=[]
		for jj in n.arange(len(eboss[around_eboss])):
			if n_cp[jj] > 0 :
				all_ages = n.array([ eboss['Chabrier_MILES_age_ssp_'+str(ii)][around_eboss][jj] for ii in n.arange(n_cp[jj]) ])
				all_masses = n.array([ eboss['Chabrier_MILES_stellar_mass_ssp_'+str(ii)][around_eboss][jj] for ii in n.arange(n_cp[jj]) ])
				sfr_inst = all_masses / all_ages
				youngest_ages.append(n.min(all_ages))
				highest_sfrs.append(n.max(sfr_inst))
			
		highest_sfrs = n.array(highest_sfrs)
		youngest_ages = n.array(youngest_ages)
		age_eboss = youngest_ages# highest_sfrs

		#age_eboss = eboss['Chabrier_MILES_stellar_mass'][around_eboss]/eboss['Chabrier_MILES_age_lightW'][around_eboss]
		# age_eboss = eboss['Chabrier_MILES_age_lightW'][around_eboss]

		position_eboss = coord.ICRS(ra_eboss*u.degree, dec_eboss*u.degree)
		sep_sky_r200c_mpc = (center.separation(position_eboss)/(cc['R200C_DEG']*u.degree)).value * r200_mpc
		sep_los_r200c_mpc = abs(aa.comoving_distance(z_eboss)-aa.comoving_distance(cc['SCREEN_CLUZSPEC'])).value
		eboss_sep_3d_per_r200 = (sep_sky_r200c_mpc**2.+sep_los_r200c_mpc**2.)**(0.5)/r200_mpc

		p.plot(eboss_sep_3d_per_r200, age_eboss, 'k,', rasterized=True)

p.xlabel('r/r200c')
#p.ylabel('stellar age [yr]')
p.ylabel('youngest age [yr]')
p.xscale('log')
p.yscale('log')
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'r200-youngestAge-z005-z01.png'))
#p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'r200-age.png'))
p.clf()

sys.exit()

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

