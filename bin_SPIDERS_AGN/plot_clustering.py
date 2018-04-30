import matplotlib.pyplot as p
import numpy as n
import os
 
#repo = os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target')
#lrg = n.loadtxt('/data36s/comparat/SDSS/LSS/LRG/lowz_boss/galaxy.2pcf', unpack=True)     
#qso = n.loadtxt('/data36s/comparat/SDSS/LSS/BigMD_QSO/BigMDPL-QSOZ_Y1Q_1481.43deg2.2pcf', unpack=True)
#cmass = n.loadtxt('/data36s/comparat/SDSS/LSS/LRG/cmass_boss/DR12v5_CMASS_North_rd0.monopole.2pcf', unpack=True)

p.figure(0, (6,6))
 
path = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_C1/C1_eRo_AGN_0.0_0.3656708956205386.2pcf'
s430 = n.loadtxt(path , unpack=True)
#p.fill_between(s430 [0], y1=s430 [1] - s430 [1]*(s430[2])**(-0.5),  y2=s430 [1] + s430 [1]*(s430[2])**(-0.5), label='eRoAGN z<0.36', color='m',alpha=0.5)
p.errorbar(s430[0], s430[1], yerr=s430[1]*(s430[2])**(-0.5), xerr=2, label='mock eRoAGN AGN')

path = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_C1/C1_RASS_AGN_0.0_0.3656708956205386.2pcf'
s430 = n.loadtxt(path , unpack=True)
#p.fill_between(s430 [0], y1=s430 [1] - s430 [1]*(s430[2])**(-0.5),  y2=s430 [1] + s430 [1]*(s430[2])**(-0.5), label='RASS AGN z<0.36', color='b',alpha=0.5)

p.errorbar(s430[0], s430[1], yerr=s430[1]*(s430[2])**(-0.5), xerr=2, label=r'mock RASS AGN ($FX>10^{-12.5}$)')

path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/clustering_agn_RL_0.005.2pcf"
d440_0 = n.loadtxt(path , unpack=True)
p.errorbar(d440_0[0], d440_0[1], yerr=d440_0[1]*(d440_0[2])**(-0.5), xerr=2, label='SPIDERS RASS')

p.title(' z<0.36')
p.xlabel('s Mpc/h')
p.ylabel('xi(s)')
#p.xscale('log')
p.yscale('log')
p.ylim((0.02,20))
p.xlim((0.,30))
p.legend(frameon=False, loc=0)
p.grid()
p.savefig(os.path.join(os.environ['HOME'], 'wwwDir', 'eRoMok/clustering/data', 'clustering_FX_lim_2018.png'))
p.clf()


p.figure(0, (6,6))
 
path = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_C2/C2_eRo_AGN_0.32342277247492107_0.7739457644425549.2pcf'
s430 = n.loadtxt(path , unpack=True)
#p.fill_between(s430 [0], y1=s430 [1] - s430 [1]*(s430[2])**(-0.5),  y2=s430 [1] + s430 [1]*(s430[2])**(-0.5), label='eRoAGN z<0.36', color='m',alpha=0.5)
p.errorbar(s430[0], s430[1], yerr=s430[1]*(s430[2])**(-0.5), label='mock eRoAGN AGN')

path = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_C2/C2_RASS_AGN_0.32342277247492107_0.7739457644425549.2pcf'
s430 = n.loadtxt(path , unpack=True)
#p.fill_between(s430 [0], y1=s430 [1] - s430 [1]*(s430[2])**(-0.5),  y2=s430 [1] + s430 [1]*(s430[2])**(-0.5), label='RASS AGN z<0.36', color='b',alpha=0.5)

p.errorbar(s430[0], s430[1], yerr=s430[1]*(s430[2])**(-0.5), label='mock RASS AGN')

path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/clustering_agn_RL_0.36_0.77_0.005.2pcf"
d440_0 = n.loadtxt(path , unpack=True)
p.errorbar(d440_0[0], d440_0[1], yerr=d440_0[1]*(d440_0[2])**(-0.5), label='data, flux limit RASS')

p.title('0.35<z<0.77')
p.xlabel('s Mpc/h')
p.ylabel('xi(s)')
p.xscale('log')
p.yscale('log')
p.ylim((0.001,10))
p.xlim((0.5,120))
p.legend(frameon=False, loc=0)
p.grid()
p.savefig(os.path.join(os.environ['HOME'], 'wwwDir', 'eRoMok/clustering/data', 'clustering_FX_lim_2018_C2.png'))
p.clf()


p.figure(0, (6,6))
 
path = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_C3/C3_eRo_AGN_0.7417971791184751_1.3449155875649825.2pcf'
s430 = n.loadtxt(path , unpack=True)
#p.fill_between(s430 [0], y1=s430 [1] - s430 [1]*(s430[2])**(-0.5),  y2=s430 [1] + s430 [1]*(s430[2])**(-0.5), label='eRoAGN z<0.36', color='m',alpha=0.5)
p.errorbar(s430[0], s430[1], yerr=s430[1]*(s430[2])**(-0.5), label='eRoAGN AGN')

path = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_C3/C3_RASS_AGN_0.7417971791184751_1.3449155875649825.2pcf'
s430 = n.loadtxt(path , unpack=True)
#p.fill_between(s430 [0], y1=s430 [1] - s430 [1]*(s430[2])**(-0.5),  y2=s430 [1] + s430 [1]*(s430[2])**(-0.5), label='RASS AGN z<0.36', color='b',alpha=0.5)

p.errorbar(s430[0], s430[1], yerr=s430[1]*(s430[2])**(-0.5), label='RASS AGN')

path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/clustering_agn_RL_0.77_1.34_0.005.2pcf"
d440_0 = n.loadtxt(path , unpack=True)
p.errorbar(d440_0[0], d440_0[1], yerr=d440_0[1]*(d440_0[2])**(-0.5), label='only flux limit RASS')

p.title('0.77<z<1.35')
p.xlabel('s Mpc/h')
p.ylabel('xi(s)')
p.xscale('log')
p.yscale('log')
p.ylim((0.001,10))
p.xlim((0.5,120))
p.legend(frameon=False, loc=0)
p.grid()
p.savefig(os.path.join(os.environ['HOME'], 'wwwDir', 'eRoMok/clustering/data', 'clustering_FX_lim_2018_C3.png'))
p.clf()

