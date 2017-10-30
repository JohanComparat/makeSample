import matplotlib.pyplot as p
import numpy as n
import os

repo = os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target')

lrg = n.loadtxt('/data36s/comparat/SDSS/LSS/LRG/lowz_boss/galaxy.2pcf', unpack=True)     
qso = n.loadtxt('/data36s/comparat/SDSS/LSS/BigMD_QSO/BigMDPL-QSOZ_Y1Q_1481.43deg2.2pcf', unpack=True)
cmass = n.loadtxt('/data36s/comparat/SDSS/LSS/LRG/cmass_boss/DR12v5_CMASS_North_rd0.monopole.2pcf', unpack=True)

path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/2rxs_all/clustering_agn_RL_0.005_LX_440.2pcf"
d440_0 = n.loadtxt(path , unpack=True) 
path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/2rxs_EXI10/clustering_agn_RL_0.005_LX_440.2pcf"
d440_1 = n.loadtxt(path , unpack=True)
path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/2rxs_EXI10_pany05/clustering_agn_RL_0.005_LX_440.2pcf"
d440_2 = n.loadtxt(path , unpack=True)

path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/2rxs_all/clustering_agn_RL_0.005_LX_443.2pcf"
d443_0 = n.loadtxt(path , unpack=True)  
path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/2rxs_EXI10/clustering_agn_RL_0.005_LX_443.2pcf"
d443_1 = n.loadtxt(path , unpack=True)
path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/2rxs_EXI10_pany05/clustering_agn_RL_0.005_LX_443.2pcf"
d443_2 = n.loadtxt(path , unpack=True)

path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/2rxs_all/clustering_agn_RL_0.005_LX_446.2pcf"
d446_0 = n.loadtxt(path , unpack=True)  
path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/2rxs_EXI10/clustering_agn_RL_0.005_LX_446.2pcf"
d446_1 = n.loadtxt(path , unpack=True)
path = "/data36s/comparat/SDSS/dr14/spiders/clustering_measurements/2rxs_EXI10_pany05/clustering_agn_RL_0.005_LX_446.2pcf"
d446_2 = n.loadtxt(path , unpack=True)
 
ref = n.loadtxt(os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/clustering_measurements/elg_favole_2016.monopole'), unpack=True)

gf440=0.89091352990878625
gf443=0.84990345086368468
gf446=0.80704933209279583


p.figure(0, (6,6))
p.errorbar(d440_0 [0], d440_0 [1], yerr=d440_0 [1]*(d440_0[2])**(-0.5), label='LX>44.0, z=0.22, all')
p.errorbar(d440_1 [0], d440_1 [1], yerr=d440_1 [1]*(d440_1[2])**(-0.5), label='LX>44.0, z=0.22, exi>10')
p.errorbar(d440_2 [0], d440_2 [1], yerr=d440_2 [1]*(d440_2[2])**(-0.5), label='LX>44.0, z=0.22, p_any>0.5')
p.plot(lrg [0], lrg [1], label='LRG z=0.28')
p.plot(cmass [0], cmass [1], label='CMASS z=0.54')
p.xlabel('s Mpc/h')
p.ylabel('xi(s)')
p.xscale('log')
p.yscale('log')
p.ylim((0.01,100))
p.xlim((0.8,100))
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(os.environ['HOME'], 'wwwDir', 'eRoMok/clustering/data', 'clustering_LX440.png'))
p.clf()

p.figure(0, (6,6))
p.errorbar(d443_0 [0], d443_0 [1], yerr=d443_0 [1]*(d443_0[2])**(-0.5), label='LX>44.3, z=0.31, all')
p.errorbar(d443_1 [0], d443_1 [1], yerr=d443_1 [1]*(d443_1[2])**(-0.5), label='LX>44.3, z=0.30, exi>10')
p.errorbar(d443_2 [0], d443_2 [1], yerr=d443_2 [1]*(d443_2[2])**(-0.5), label='LX>44.3, z=0.30, p_any>0.5')
p.plot(lrg [0], lrg [1], label='LRG z=0.28')
p.plot(cmass [0], cmass [1], label='CMASS z=0.54')
p.xlabel('s Mpc/h')
p.ylabel('xi(s)')
p.xscale('log')
p.yscale('log')
p.ylim((0.01,100))
p.xlim((0.8,100))
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(os.environ['HOME'], 'wwwDir', 'eRoMok/clustering/data', 'clustering_LX443.png'))
p.clf()


p.figure(0, (6,6))
p.errorbar(d446_0 [0], d446_0 [1], yerr=d446_0 [1]*(d446_0[2])**(-0.5), label='LX>44.6, z=0.41, all')
p.errorbar(d446_1 [0], d446_1 [1], yerr=d446_1 [1]*(d446_1[2])**(-0.5), label='LX>44.6, z=0.40, exi>10')
p.errorbar(d446_2 [0], d446_2 [1], yerr=d446_2 [1]*(d446_2[2])**(-0.5), label='LX>44.6, z=0.40, p_any>0.5')
p.plot(lrg [0], lrg [1], label='LRG z=0.28')
p.plot(cmass [0], cmass [1], label='CMASS z=0.54')
p.xlabel('s Mpc/h')
p.ylabel('xi(s)')
p.xscale('log')
p.yscale('log')
p.ylim((0.01,100))
p.xlim((0.8,100))
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(os.environ['HOME'], 'wwwDir', 'eRoMok/clustering/data', 'clustering_LX446.png'))
p.clf()



