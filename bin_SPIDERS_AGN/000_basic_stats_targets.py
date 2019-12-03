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
# AGN
path_2_targets = os.path.join(target_dir, 'spiderstargetSequelsAGN-SPIDERS_RASS_AGN-v1.1.fits'     ) 
tgA_RASS_s = Selection()
tgA_RASS_s.tg = fits.open(path_2_targets)#[1].data

path_2_targets = os.path.join(target_dir, 'spiderstargetAGN-SPIDERS_XMMSL_AGN-v3.1.fits'           )   
tgA_XMMSL = Selection()
tgA_XMMSL.tg = fits.open(path_2_targets)#[1].data

path_2_targets = os.path.join(target_dir, 'spiderstargetAGN-SPIDERS_RASS_AGN-v2.1.fits'            )   
tgA_RASS = Selection()
tgA_RASS.tg = fits.open(path_2_targets)#[1].data

# CLUSTER galaxies
path_2_targets = os.path.join(target_dir, 'spiderstargetSequelsClus-SPIDERS_RASS_CLUS-v1.0.fits'   )   
tgC_RASS_s = Selection()
tgC_RASS_s.tg = fits.open(path_2_targets)#[1].data

path_2_targets = os.path.join(target_dir, 'spiderstargetClusters-SPIDERS_XCLASS_CLUS-v1.1.fits'    )   
tgC_XCLASS = Selection()
tgC_XCLASS.tg = fits.open(path_2_targets)#[1].data

path_2_targets = os.path.join(target_dir, 'spiderstargetClusters-SPIDERS_RASS_CLUS-v1.1.fits'      )   
tgC_RASS = Selection()
tgC_RASS.tg = fits.open(path_2_targets)#[1].data

# 2RXS
path_2_2RXS_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits') # '2RXS_AllWISE_catalog_paper_2017May26_v5_11_0_sdss_26_VERON_MASKED_GAIA_star_mask.fits')
full_2RXS = fits.open(path_2_2RXS_cat)[1].data

path_2_9028_spec = os.path.join(catalog_dir, 'spiderstargetAGN-SPIDERS_RASS_AGN-v2.1_v5_13_0_sdss_26_VERON_2QZ_MASKED.fits')
full_9028_spec = fits.open(path_2_9028_spec)[1].data

path_2_819_spec = os.path.join(catalog_dir, 'spiderstargetAGN-SPIDERS_XMMSL_AGN-v3.1_v5_13_0_sdss_26_VERON_2QZ_MASKED.fits')
full_819_spec = fits.open(path_2_819_spec)[1].data

# catalogs available :
#path_2_cat = os.path.join(catalog_dir, 'results_all.fits.gz')
#path_2_cat = os.path.join(catalog_dir, 'results_VAC_2RXS.fits')
path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_2RXS_DR16.fits')
VAC_2RXS = fits.open(path_2_cat)[1].data

path_2_cat = os.path.join(catalog_dir, 'VAC_SPIDERS_XMMSL2_DR16.fits')
VAC_XMMSL = fits.open(path_2_cat)[1].data

path_2_cat = os.path.join(catalog_dir, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits')
tg_VAC_all_2RXS = Selection()
tg_VAC_all_2RXS.tg = fits.open(path_2_cat) #fits.open(path_2_targets)#[1].data

path_2_cat = os.path.join(catalog_dir, 'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI.fits')
tg_VAC_all_XMMSL = Selection()
tg_VAC_all_XMMSL.tg = fits.open(path_2_cat) #fits.open(path_2_targets)#[1].data


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
#masks[os.path.basename(mangle_DR7_file)]   = pymangle.Mangle(mangle_DR7_file )
masks[os.path.basename(mangleSPIDERS_dr16S_file )]  = pymangle.Mangle(mangleSPIDERS_dr16S_file  )

mask_labels = {}
mask_labels[os.path.basename(mangle_DR16_file )] = 'DR16'
mask_labels[os.path.basename(mangle_DR12N_file)] = 'DR12'
mask_labels[os.path.basename(mangle_DR12S_file)] = 'DR12S'
mask_labels[os.path.basename(mangle_DR14N_file)] = 'DR14N'
mask_labels[os.path.basename(mangle_DR14S_file)] = 'DR14S'
#mask_labels[os.path.basename(mangle_DR7_file)] ='DR7'
mask_labels[os.path.basename(mangleSPIDERS_dr16S_file )] = 'DR16S'


def get_in_mask(masks= masks, targets = tgA_RASS_s.tg):
	in_masks = {}
	for kk in masks.keys():
		in_masks[kk] = masks[kk].contains(targets[2].data['RA'], targets[2].data['DEC'])  
	return in_masks

def get_in_mask_2(masks= masks, targets = tgA_RASS_s.tg):
	in_masks = {}
	for kk in masks.keys():
		in_masks[kk] = masks[kk].contains(targets[1].data['ALLW_RA'], targets[1].data['ALLW_DEC'])  
	return in_masks


def get_N_in_mask(in_masks):
	N_in_masks = {}
	for kk in in_masks.keys():
		N_in_masks[kk] = len(in_masks[kk].nonzero()[0]) 
	return N_in_masks

#masks[os.path.basename(mangle_DR7_file )].get_polyids( tgA_RASS_s.tg[2].data['RA'], tgA_RASS_s.tg[2].data['DEC'] )


tg_VAC_all_XMMSL.in_masks = get_in_mask_2(masks, tg_VAC_all_XMMSL.tg)
tg_VAC_all_XMMSL.N = tg_VAC_all_XMMSL.tg[1].header['NAXIS2']
tg_VAC_all_XMMSL.N_in_masks = get_N_in_mask(tg_VAC_all_XMMSL.in_masks)
print('tg_VAC_all_XMMSL', tg_VAC_all_XMMSL.N, tg_VAC_all_XMMSL.N_in_masks)

tg_VAC_all_2RXS.in_masks = get_in_mask_2(masks, tg_VAC_all_2RXS.tg)
tg_VAC_all_2RXS.N = tg_VAC_all_2RXS.tg[1].header['NAXIS2']
tg_VAC_all_2RXS.N_in_masks = get_N_in_mask(tg_VAC_all_2RXS.in_masks)
print('tg_VAC_all_2RXS', tg_VAC_all_2RXS.N, tg_VAC_all_2RXS.N_in_masks)

tgA_RASS_s.in_masks = get_in_mask(masks, tgA_RASS_s.tg)
tgA_RASS_s.N = tgA_RASS_s.tg[2].header['NAXIS2']
tgA_RASS_s.N_in_masks = get_N_in_mask(tgA_RASS_s.in_masks)
print('tgA_RASS_s', tgA_RASS_s.N, tgA_RASS_s.N_in_masks)

tgA_XMMSL.in_masks = get_in_mask(masks, tgA_XMMSL.tg)
tgA_XMMSL.N = tgA_XMMSL.tg[2].header['NAXIS2']
tgA_XMMSL.N_in_masks = get_N_in_mask(tgA_XMMSL.in_masks)
print('tgA_XMMSL', tgA_XMMSL.N, tgA_XMMSL.N_in_masks)

tgA_RASS.in_masks = get_in_mask(masks, tgA_RASS.tg)
tgA_RASS.N = tgA_RASS.tg[2].header['NAXIS2']
tgA_RASS.N_in_masks = get_N_in_mask(tgA_RASS.in_masks)
print('tgA_RASS', tgA_RASS.N, tgA_RASS.N_in_masks)

tgC_RASS_s.in_masks = get_in_mask(masks, tgC_RASS_s.tg)
tgC_RASS_s.N = tgC_RASS_s.tg[2].header['NAXIS2']
tgC_RASS_s.N_in_masks = get_N_in_mask(tgC_RASS_s.in_masks)
print('tgC_RASS_s', tgC_RASS_s.N, tgC_RASS_s.N_in_masks)

tgC_RASS.in_masks = get_in_mask(masks, tgC_RASS.tg)
tgC_RASS.N = tgC_RASS.tg[2].header['NAXIS2']
tgC_RASS.N_in_masks = get_N_in_mask(tgC_RASS.in_masks)
print('tgC_RASS', tgC_RASS.N, tgC_RASS.N_in_masks)

tgC_XCLASS.in_masks = get_in_mask(masks, tgC_XCLASS.tg)
tgC_XCLASS.N = tgC_XCLASS.tg[2].header['NAXIS2']
tgC_XCLASS.N_in_masks = get_N_in_mask(tgC_XCLASS.in_masks)
print('tgC_XCLASS', tgC_XCLASS.N, tgC_XCLASS.N_in_masks)


# determine the number with spectra in 5_13_0 
kk = 'DR16_wSEQUELS_footprint_clipped_complete.pol'

wspec_819 = masks[kk].contains(full_819_spec['ALLWISE_RA'][full_819_spec['in_BOSS_v5_13_0'] ], full_819_spec['ALLWISE_DEC'][full_819_spec['in_BOSS_v5_13_0'] ])
wspec_9028 = masks[kk].contains(full_9028_spec['ALLWISE_RA'][full_9028_spec['in_BOSS_v5_13_0'] ], full_9028_spec['ALLWISE_DEC'][full_9028_spec['in_BOSS_v5_13_0'] ])

tgA_XMMSL.Nspec = len(wspec_819.nonzero()[0])
tgA_RASS.Nspec = len(wspec_9028.nonzero()[0])

print( tgA_RASS.Nspec,tgA_RASS.N_in_masks[kk]  ,tgA_RASS.Nspec/tgA_RASS.N_in_masks[kk]   ) 
print( tgA_XMMSL.Nspec,tgA_XMMSL.N_in_masks[kk],tgA_XMMSL.Nspec/tgA_XMMSL.N_in_masks[kk] )

RRs = {}
#for kk in masks.keys():
kk = 'mask_DR12v5_CMASS_North.ply'
RRs[kk] = masks[kk].genrand(20000)  
kk = 'DR16_wSEQUELS_footprint_clipped_complete.pol'
RRs[kk] = masks[kk].genrand(10000)  

fig_out = os.path.join(figure_dir, 'ra-dec-mask.png' )

p.figure(1, (5.5,5.5))
p.axes([0.15, 0.15, 0.8, 0.77])
p.tight_layout()
kk = 'mask_DR12v5_CMASS_North.ply'
p.plot(RRs[kk][0], RRs[kk][1], 'ko', ls='', rasterized=True )
kk = 'DR16_wSEQUELS_footprint_clipped_complete.pol'
p.plot(RRs[kk][0], RRs[kk][1], 'ro', ls='', label='eBOSS area', rasterized=True )
p.xlabel('R.A. [deg]')
p.ylabel('Declination [deg]')
p.grid()
#p.legend(frameon=True, loc=0)
p.savefig(fig_out)
p.clf()

#kk = 'mask_DR12v5_CMASS_North.ply'
#rxs_indr13 = masks[kk].contains(full_2RXS['ALLW_RA'],full_2RXS['ALLW_DEC'])
#print( len(rxs_indr13.nonzero()[0])  )

#kk = 'DR16_wSEQUELS_footprint_clipped_complete.pol'
#indr16 = masks[kk].contains(full_2RXS['ALLW_RA'],full_2RXS['ALLW_DEC'])
#sel = full_2RXS['NWAY_p_any']>=0.01 
#rxs_indr16 = (sel) & (indr16)
#print( len(indr16.nonzero()[0]), len(rxs_indr16.nonzero()[0]) )