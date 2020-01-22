"""
python plot_XLF.py

"""
from catalog_lib import *

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

env = "MD10" # sys.argv[1]

area_eboss_dr16 = 5128.

sig_fx = 0.1
MAG_MIN_LIM = 15
MAG_LIM = 22.5
EXI_ML_min = 6.5

# LF binning
dlogf = 0.25
fbins = np.arange(36,47,dlogf)
xf = np.arange(36,47,dlogf)[:-1]+dlogf/2.

redshift = data_RXS.data['Z_BEST']
redshift[redshift<0]=0

dz=0.2
z_min=0.0 
z_max = np.max(redshift[agnZ]) # np.max([np.max(redshift),np.max(Z_XMMSL2)])
zs= np.arange(z_min, z_max + dz, dz)
xzs=(zs[1:]+zs[:-1])/2.

dataDir = os.path.join(os.environ['GIT_VS'], 'data')
xlf_dir = os.path.join(dataDir, 'LF_SMF/LXFunction/aird_2015/')
pta = lambda name : os.path.join(xlf_dir, name)
pta('XLF_soft_00_02.csv')
pta('XLF_soft_04_06.csv')
pta('XLF_soft_08_10.csv')
pta('XLF_soft_15_20.csv')
pta('XLF_soft_25_35.csv')
x_aj12, y_aj12 = np.loadtxt(os.path.join(dataDir, 'LF_SMF', 'LXFunction', 'ajello_12.csv'), delimiter=',', unpack=True)

log_LX_hard, log_dndv_up,  log_dndv_low = np.loadtxt(os.path.join(dataDir, 'LF_SMF', 'LXFunction', 'buchner_2015.txt'), unpack=True)

z0, z1, Lxmin, Lxmax, Nobj, Nmdl, phi_a = np.loadtxt(os.path.join(os.environ['GIT_VS'], 'data/LF_SMF', 'LXFunction', 'hasinger_2005_soft_XLF.ascii'), unpack=True)
z_center = ( z0 + z1 ) * 0.5
Lx_c = ( Lxmin + Lxmax ) * 0.5 - 2*np.log10(0.7)
phi = phi_a * 0.7**3.

# aird model

# hard
kz_h = lambda z : 10**(-4.03 - 0.19*(1+z))
Ls_h = lambda z : 10**(44.84 - np.log10( ((1+2.0)/(1+z))**3.87 + ((1+2.0)/(1+z))**(-2.12) ) )
phi_h = lambda L, z  : kz_h(z) / ( (L/Ls_h(z))**0.48 + (L/Ls_h(z))**2.27)

# soft
kz_s = lambda z : 10**(-4.28 - 0.22*(1+z))
Ls_s = lambda z : 10**(44.93 - np.log10( ((1+2.31)/(1+z))**3.39 + ((1+2.31)/(1+z))**(-3.58) ) )
phi_s = lambda L, z  : kz_s(z) / ( (L/Ls_s(z))**0.44 + (L/Ls_s(z))**2.18)


def get_nzs(lx, z, fbins = fbins, zmin=0, zmax=6, area=area_eboss_dr16):
  vol = area * (cosmoMD.comoving_volume(zmax).value-cosmoMD.comoving_volume(zmin).value) * np.pi/129600.
  zsel = (z>=zmin)&(z<zmax)&(lx>0)
  nall = np.histogram(lx[zsel]   , fbins)[0]/vol/dlogf/np.log(10)
  nallN = np.histogram(lx[zsel]  , fbins)[0]
  return nall, nallN, vol

# AGN NZ from the mock catalogue
path_2_mock = os.path.join(os.environ[env],'SPIDERS_AGN_all.fit')

print(path_2_mock)
hd_mock_a = fits.open(path_2_mock)[1].data
fx_true_a = hd_mock_a['agn_FX_soft'] 
simulation_area = 30575. #129600/np.pi
hd_mock = hd_mock_a[ ( fx_true_a > 2e-15 ) & ( abs(hd_mock_a['g_lat']) > 15 ) ]
fx_true = hd_mock['agn_FX_soft'] 
redshift_R = hd_mock['redshift_R']
LX_soft = hd_mock['AGN_LX_soft']
N_obj = len(fx_true)

for zmin in zs[:-1]:
	print(zmin, dz)
	zmax = zmin+dz
	z_mean = (zmin+zmax)*0.5*np.ones_like(xf)
	ok_hasinger = (z_center>zmin)&(z_center<=zmax)
	
	p.figure(1, (6,6))
	p.axes([0.2,0.18,0.75,0.75])
	p.plot(xf,phi_h(10**xf,z_mean), c='cyan' , ls='dashed', lw=2, label='Ai15')#Aird 2-10 keV LADE')
	#p.fill_between(xf, y1=(nall)*(1-(nallN)**(-0.5)), y2=(nall)*(1+(nallN)**(-0.5)), color='r', alpha=0.7, label='Soft mock')# 0.5-2 keV')
	if len(Lx_c[ok_hasinger])>0:
		p.fill_between(Lx_c[ok_hasinger], y1=phi[ok_hasinger]*(1-(Nobj[ok_hasinger])**(-0.5)), y2=phi[ok_hasinger]*(1+(Nobj[ok_hasinger])**(-0.5)), color='G', alpha=0.7, label='Ha05', lw=2)
	if z_mean[0]>1 and z_mean[0]<1.5:
		p.fill_between(log_LX_hard, y1=10**log_dndv_low, y2=10**log_dndv_up, color='g', alpha=0.7, label='Bu15', lw=2)

	log_LX_hard, log_dndv_up,  log_dndv_low
	#nall, nallN, vol = get_nzs(LX_soft, redshift_R, fbins, zmin, zmax, simulation_area)
	#p.fill_between(xf,  y1=(nall)*(1-(nallN)**(-0.5)), y2=(nall)*(1+(nallN)**(-0.5)), color='g', alpha=0.7, label='Mock', lw=2)# 2-10  keV')

	nall, nallN, vol = get_nzs(data_RXS.data['RXS_LOGLX'][data_RXS.observed], data_RXS.data['Z_BEST'][data_RXS.observed], fbins, zmin, zmax, area_eboss_dr16)
	p.fill_between(xf,  y1=(nall)*(1-(nallN)**(-0.5)), y2=(nall)*(1+(nallN)**(-0.5)), color='r', alpha=0.7, label='2RXS', lw=2)# 2-10  keV')

	nall, nallN, vol = get_nzs(data_XMM.data['XMMSL2_LOGLX_FULL'][data_XMM.observed], data_XMM.data['Z_BEST'][data_XMM.observed], fbins, zmin, zmax, area_eboss_dr16)
	p.fill_between(xf,  y1=(nall)*(1-(nallN)**(-0.5)), y2=(nall)*(1+(nallN)**(-0.5)), color='b', alpha=0.7, label='XMMSL2', lw=2)# 2-10  keV')

	p.xlabel(r'$\log_{10}(L_X/[erg/s])$')
	p.ylabel(r'$\Phi$=dN/dlogL/dV [1/Mpc$^3$/dex]')
	#if zmin<0.1:
	p.legend(loc=3, fontsize=14)
	p.yscale('log')
	p.xlim((42., 47.5))
	p.ylim((1/(2*vol),1e-3))
	p.title(str(np.round(zmin,2))+"<z<"+str(np.round(zmax,2)))
	p.grid()
	p.savefig(os.path.join(figure_dir, "XLF_soft_"+str(np.round(zmin,2))+"_"+str(np.round(zmax,2))+".png"))
	p.clf()
