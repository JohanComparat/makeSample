import numpy as np
import sys,os


def write_cats(z_mask, name):
	path_2_elg_cat = os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'ebosselg.JC28Nov2017_data' )
	ra, dec, z, w, plate, mjd, fiberid, sector, tsr, pssr = np.loadtxt(path_2_elg_cat, unpack=True)
	z_int = (z*10000).astype('int')
	keep_in_catalog = (np.in1d(z_int, z_mask)==False)&(dec<5)
	np.savetxt(path_2_elg_cat+name, np.transpose([ra[keep_in_catalog], dec[keep_in_catalog], z[keep_in_catalog], w[keep_in_catalog] ]))
	path_2_elg_cat = os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'ebosselg.JC28Nov2017_rand' )
	ra, dec, z, w, sector, tsr, pssr = np.loadtxt(path_2_elg_cat, unpack=True)
	z_int = (z*10000).astype('int')
	keep_in_catalog = (np.in1d(z_int, z_mask)==False)&(dec<5)
	np.savetxt(path_2_elg_cat+name, np.transpose([ra[keep_in_catalog], dec[keep_in_catalog], z[keep_in_catalog], w[keep_in_catalog] ]))

z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2O3Hb_th_3.txt'))*10000).astype('int')
write_cats(z_mask, name='.O2O3Hb_th_3.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2O3_th_3.txt'))  *10000).astype('int')
write_cats(z_mask, name='.O2O3_th_3.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2_th_3.txt'))    *10000).astype('int')
write_cats(z_mask, name='.O2_th_3.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2O3Hb_th_2.txt'))*10000).astype('int')
write_cats(z_mask, name='.mask_O2O3Hb_th_2.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2O3_th_2.txt'))  *10000).astype('int')
write_cats(z_mask, name='.mask_O2O3_th_2.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2_th_2.txt'))    *10000).astype('int')
write_cats(z_mask, name='.mask_O2_th_2.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2O3Hb_th_1.txt'))*10000).astype('int')
write_cats(z_mask, name='.mask_O2O3Hb_th_1.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2O3_th_1.txt'))  *10000).astype('int')
write_cats(z_mask, name='.mask_O2O3_th_1.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2_th_1.txt'))    *10000).astype('int')
write_cats(z_mask, name='.mask_O2_th_1.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2O3Hb_th_1.5.txt'))*10000).astype('int')
write_cats(z_mask, name='.mask_O2O3Hb_th_1.5.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2O3_th_1.5.txt'))  *10000).astype('int')
write_cats(z_mask, name='.mask_O2O3_th_1.5.txt')
z_mask = ( np.loadtxt(os.path.join(os.environ['HOME'], 'data/SDSS/eBOSS/ELG', 'elg_z_mask_O2_th_1.5.txt'))    *10000).astype('int')
write_cats(z_mask, name='.mask_O2_th_1.5.txt')

# then write the catalog with the 
