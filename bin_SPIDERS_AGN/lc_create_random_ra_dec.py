import numpy as n

# about 20e6 full sky randoms :
size=2000000
for i in range(20):
	uu = n.random.uniform(size=size)
	dec = n.arccos(1-2*uu)* 180/n.pi - 90.
	ra = n.random.uniform(size=size)*2*n.pi* 180/n.pi
	topdir = '/data44s/eroAGN_WG_DATA/DATA/randoms/'
	n.savetxt(topdir+'random-ra-dec-'+str(i)+'.txt', n.transpose([ra,dec]), header="RA DEC")
