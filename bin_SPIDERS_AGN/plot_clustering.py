import matplotlib.pyplot as p
import numpy as n
import os

repo = os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target')
     

d_N_005 = n.loadtxt(os.path.join(repo,'clustering_agn_N_RL_0.005.2pcf' ), unpack=True)  
d_N_0216 = n.loadtxt(os.path.join(repo,'clustering_agn_N_RL_0.0216.2pcf'), unpack=True)  
d_S_005 = n.loadtxt(os.path.join(repo,'clustering_agn_S_RL_0.005.2pcf' ), unpack=True)  
d_S_015 = n.loadtxt(os.path.join(repo,'clustering_agn_S_RL_0.015.2pcf' ), unpack=True)

ref = n.loadtxt(os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/clustering/elg_favole_2016.monopole'), unpack=True)

p.figure(0, (6,6))
#p.axes()
p.plot(d_N_005 [0], d_N_005 [1], label='N_005 ')
p.plot(d_N_0216[0], d_N_0216[1], label='N_0216')
p.plot(d_S_005 [0], d_S_005 [1], label='S_005 ')
p.plot(d_S_015 [0], d_S_015 [1], label='S_015 ')
p.plot(ref [0], ref [1], label='ELG ')

p.xlabel('s Mpc/h')
p.ylabel('xi(s)')
p.xscale('log')
p.yscale('log')
p.xlim((0.1,100))
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(repo, 'clustering.png'))
p.clf()
