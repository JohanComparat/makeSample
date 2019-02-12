import numpy as n

# about 20e6 full sky randoms :
size=2000000
for i in range(20):
	x = (n.random.uniform(size=size)-0.5)*4.#*u_max
	y = (n.random.uniform(size=size)-0.5)*4.#*v_max
	z = (n.random.uniform(size=size)-0.5)*4.#*u_max
	r2 = x*x + y*y + z*z
	sel = (r2>1.)&(r2<=2)
	dec  = n.arccos( z / r2 ) * 180/n.pi - 90.
	ra = n.arctan2( y , x ) * 180/n.pi + 180.
	topdir = '/data44s/eroAGN_WG_DATA/DATA/randoms/'
	n.savetxt(topdir+'random-ra-dec-'+str(i)'.txt', n.transpose([ra[sel],dec[sel]]), header="RA DEC")
