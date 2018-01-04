import numpy as np
import sys,os

maskLambda = np.loadtxt(os.path.join(os.environ['GIT_FF'],'data',"dr12-sky-mask.txt"), unpack=True)

z_array = np.arange(0.6,1.2,0.0001)
o2a_wl = 3726.032 * (1+z_array)
o2b_wl = 3728.814 * (1+z_array)
o3b_wl = 5006.841 * (1+z_array)
hb_wl  = 4861.331 * (1+z_array)

def get_mask(wl,margin = 1.5):
	ratio = np.min(abs(10000.*np.log10(np.outer(wl, 1./maskLambda))), axis=1)
	return ratio <= margin


mask_o2a = get_mask(o2a_wl, 3)
mask_o2b = get_mask(o2b_wl, 3)
mask_o3b = get_mask(o3b_wl, 3)
mask_hb  = get_mask(hb_wl , 3)

mask_all = (mask_o2a )&( mask_o2b )&( mask_o3b )&( mask_hb)
mask_o2o3 = (mask_o2a )&( mask_o2b )&( mask_o3b )
mask_o2 = (mask_o2a )&( mask_o2b )

print(len(mask_all.nonzero()[0]), len(mask_o2o3.nonzero()[0]), len(mask_o2.nonzero()[0]), len(z_array))

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2O3Hb_th_3.txt'), np.transpose([z_array[mask_all]]), fmt='%1.4f')

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2O3_th_3.txt'), np.transpose([z_array[mask_o2o3]]), fmt='%1.4f')

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2_th_3.txt'), np.transpose([z_array[mask_o2]]), fmt='%1.4f')


mask_o2a = get_mask(o2a_wl, 2)
mask_o2b = get_mask(o2b_wl, 2)
mask_o3b = get_mask(o3b_wl, 2)
mask_hb  = get_mask(hb_wl , 2)

mask_all = (mask_o2a )&( mask_o2b )&( mask_o3b )&( mask_hb)
mask_o2o3 = (mask_o2a )&( mask_o2b )&( mask_o3b )
mask_o2 = (mask_o2a )&( mask_o2b )

print(len(mask_all.nonzero()[0]), len(mask_o2o3.nonzero()[0]), len(mask_o2.nonzero()[0]), len(z_array))

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2O3Hb_th_2.txt'), np.transpose([z_array[mask_all]]), fmt='%1.4f')

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2O3_th_2.txt'), np.transpose([z_array[mask_o2o3]]), fmt='%1.4f')

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2_th_2.txt'), np.transpose([z_array[mask_o2]]), fmt='%1.4f')


mask_o2a = get_mask(o2a_wl, 1)
mask_o2b = get_mask(o2b_wl, 1)
mask_o3b = get_mask(o3b_wl, 1)
mask_hb  = get_mask(hb_wl , 1)

mask_all = (mask_o2a )&( mask_o2b )&( mask_o3b )&( mask_hb)
mask_o2o3 = (mask_o2a )&( mask_o2b )&( mask_o3b )
mask_o2 = (mask_o2a )&( mask_o2b )

print(len(mask_all.nonzero()[0]), len(mask_o2o3.nonzero()[0]), len(mask_o2.nonzero()[0]), len(z_array))

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2O3Hb_th_1.txt'), np.transpose([z_array[mask_all]]), fmt='%1.4f')

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2O3_th_1.txt'), np.transpose([z_array[mask_o2o3]]), fmt='%1.4f')

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2_th_1.txt'), np.transpose([z_array[mask_o2]]), fmt='%1.4f')


mask_o2a = get_mask(o2a_wl, 1.5)
mask_o2b = get_mask(o2b_wl, 1.5)
mask_o3b = get_mask(o3b_wl, 1.5)
mask_hb  = get_mask(hb_wl , 1.5)

mask_all = (mask_o2a )&( mask_o2b )&( mask_o3b )&( mask_hb)
mask_o2o3 = (mask_o2a )&( mask_o2b )&( mask_o3b )
mask_o2 = (mask_o2a )&( mask_o2b )

print(len(mask_all.nonzero()[0]), len(mask_o2o3.nonzero()[0]), len(mask_o2.nonzero()[0]), len(z_array))

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2O3Hb_th_1.5.txt'), np.transpose([z_array[mask_all]]), fmt='%1.4f')

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2O3_th_1.5.txt'), np.transpose([z_array[mask_o2o3]]), fmt='%1.4f')

np.savetxt(os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2_th_1.5.txt'), np.transpose([z_array[mask_o2]]), fmt='%1.4f')


import astropy.io.fits as fits
import numpy as np

path_2_mask = os.path.join(os.environ['HOME'], 'data', 'elg_z_mask_O2O3Hb_th_1.5.txt')
path_2_elg_cat = os.path.join(os.environ['HOME'], 'data', 'eboss21.v5_10_7.latest.fits' )

z_mask = (np.loadtxt(path_2_mask)*10000).astype('int')
hdulist = fits.open(path_2_elg_cat)

z_int = (hdulist[1].data['Z']*10000).astype('int')

keep_in_catalog = (np.in1d(z_int, z_mask)==False)

# then write the catalog with the 
new_data = hdulist[1].data[keep_in_catalog]
