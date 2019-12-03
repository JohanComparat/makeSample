"""
python2 to have access to mangle
python2 basic-stats.py

pyCONDA
python to have access to pymangle

cd ~/software/mangle2.2/bin/

pixelize -P10s /home/comparat/data/AGN_clustering/footprint/DR16_wSEQUELS_footprint_clipped_complete.ply jfhp
vim jfhp
# delete first lines
snap jfhp jfh
balkanize jfh jb
unify jb /home/comparat/data/AGN_clustering/footprint/DR16_wSEQUELS_footprint_clipped_complete.pol

pixelize -P10s /home/comparat/data/AGN_clustering/footprint/sdss_dr72safe0_res6d.pol jfhp
vim jfhp
# delete first lines
snap jfhp jfh
balkanize jfh jb
unify jb /home/comparat/data/AGN_clustering/footprint/sdss_dr72.pol

"""
from catalog_lib import *

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p


class Selection: pass

# target catalogs : 

# 2RXS on eBOSS
path_2_2RXS_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits') # '2RXS_AllWISE_catalog_paper_2017May26_v5_11_0_sdss_26_VERON_MASKED_GAIA_star_mask.fits')
full_2RXS = fits.open(path_2_2RXS_cat)[1].data

# catalogs available :
#path_2_cat = os.path.join(catalog_dir, 'results_all.fits.gz')
#path_2_cat = os.path.join(catalog_dir, 'results_VAC_2RXS.fits')
path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_2RXS_DR16.fits')
VAC_2RXS = fits.open(path_2_cat)[1].data

path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_XMMSL2_DR16.fits')
VAC_XMMSL = fits.open(path_2_cat)[1].data


masks = {}
# masks and footprint
mangle_DR12N_file = os.path.join(footprint_dir,'mask_DR12v5_CMASS_North.ply')
mangle_DR12S_file = os.path.join(footprint_dir,'mask_DR12v5_CMASS_South.ply')
mangle_DR16_file  = os.path.join(footprint_dir,'eBOSS_QSOandLRG_fullfootprintgeometry_noveto.ply')
mangle_DR14N_file = os.path.join(footprint_dir,'mask_DR14_QSO_N.ply')
mangle_DR14S_file = os.path.join(footprint_dir,'mask_DR14_QSO_S.ply')
mangle_DR7_file = os.path.join(footprint_dir,'sdss_dr72.pol')
mangleSPIDERS_dr16S_file = os.path.join(footprint_dir, 'DR16_wSEQUELS_footprint_clipped_complete.pol' )  

# read a mangle polygon file
masks[os.path.basename(mangle_DR16_file )] = pymangle.Mangle(mangle_DR16_file  ) # pymangle.Mangle
masks[os.path.basename(mangle_DR12N_file)] = pymangle.Mangle(mangle_DR12N_file )
masks[os.path.basename(mangle_DR12S_file)] = pymangle.Mangle(mangle_DR12S_file )
masks[os.path.basename(mangle_DR14N_file)] = pymangle.Mangle(mangle_DR14N_file )
masks[os.path.basename(mangle_DR14S_file)] = pymangle.Mangle(mangle_DR14S_file )
masks[os.path.basename(mangle_DR7_file)]   = pymangle.Mangle(mangle_DR7_file )
masks[os.path.basename(mangleSPIDERS_dr16S_file )]  = pymangle.Mangle(mangleSPIDERS_dr16S_file  )

mask_labels = {}
mask_labels[os.path.basename(mangle_DR16_file )] = 'DR16'
mask_labels[os.path.basename(mangle_DR12N_file)] = 'DR12'
mask_labels[os.path.basename(mangle_DR12S_file)] = 'DR12S'
mask_labels[os.path.basename(mangle_DR14N_file)] = 'DR14N'
mask_labels[os.path.basename(mangle_DR14S_file)] = 'DR14S'
mask_labels[os.path.basename(mangle_DR7_file)] ='DR7'
mask_labels[os.path.basename(mangleSPIDERS_dr16S_file )] = 'DR16S'


def get_in_mask(masks, targets):
	in_masks = {}
	for kk in masks.keys():
		in_masks[kk] = masks[kk].contains(targets[2].data['RA'], targets[2].data['DEC'])  
	return in_masks


def get_N_in_mask(in_masks):
	N_in_masks = {}
	for kk in in_masks.keys():
		N_in_masks[kk] = len(in_masks[kk].nonzero()[0]) 
	return N_in_masks


RRs = {}
#for kk in masks.keys():
kk = 'mask_DR12v5_CMASS_North.ply'
RRs[kk] = masks[kk].genrand(20000)  
kk = 'DR16_wSEQUELS_footprint_clipped_complete.pol'
RRs[kk] = masks[kk].genrand(10000)  


kk = 'mask_DR12v5_CMASS_North.ply'
rxs_indr13 = masks[kk].contains(full_2RXS['ALLW_RA'],full_2RXS['ALLW_DEC'])
print( len(rxs_indr13.nonzero()[0])  )

kk = 'DR16_wSEQUELS_footprint_clipped_complete.pol'
indr16 = masks[kk].contains(full_2RXS['ALLW_RA'],full_2RXS['ALLW_DEC'])
sel = full_2RXS['NWAY_p_any']>=0.01 
rxs_indr16 = (sel) & (indr16)
print( len(indr16.nonzero()[0]), len(rxs_indr16.nonzero()[0]) )
not_indr16 = ( indr16 == False )

fig_out = os.path.join(figure_dir, 'mask-out.png' )

p.figure(1, (5.5,5.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()
p.plot(full_2RXS['ALLW_RA'][indr16],full_2RXS['ALLW_DEC'][indr16], 'r,', ls='', rasterized=True )
p.plot(full_2RXS['ALLW_RA'][not_indr16],full_2RXS['ALLW_DEC'][not_indr16], 'k+', ls='', rasterized=True )
p.xlabel('R.A. [deg]')
p.ylabel('Declination [deg]')
p.grid()
#p.legend(frameon=True, loc=0)
p.savefig(fig_out)
p.clf()
