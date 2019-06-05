"""
python3.4 create_stack_lists_X_AGN.py
"""
import sys
import numpy as n
from scipy.interpolate import interp1d
import os
import astropy.io.fits as fits
import pickle

version = 'v1'

zmins = n.hstack((-0.01, 0.005, n.arange(0.1,3.0,0.1) ))
zmaxs = n.hstack((0.005, n.arange(0.,3.0,0.1)+0.5))
zmins = n.array([-0.005])
zmaxs = n.array([5.])

# on laptop
#path_2_cats = os.path.join(os.environ['HOME'], 'data', 'spiders', 'agn')
# on servers
path_2_cats = os.path.join(os.environ['HOME'], 'data1/SDSS/dr16/spiders')
path_2_stack_lists = os.path.join(os.environ['HOME'], 'SDSS/stacks/X_AGN')

# 4 catalogs 
path_2_XMMSL2 = os.path.join(path_2_cats, 'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR14_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')
path_2_2RXS   = os.path.join(path_2_cats, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR14_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')
path_2_XXL    = os.path.join(path_2_cats, 'Menzel16_spAll_XXL_plate.fits')
path_2_S82X   = os.path.join(path_2_cats, 'LaMassa19_spAll_S82X_plate.fits')

cat_XMMSL2 = fits.open(path_2_XMMSL2)[1].data
cat_2RXS   = fits.open(path_2_2RXS  )[1].data
cat_XXL    = fits.open(path_2_XXL   )[1].data
cat_S82X   = fits.open(path_2_S82X  )[1].data

# 

# S82X lists
survey = 'S82X'
cat = cat_S82X

z_best = cat['Z']
out = cat['S82XVI_Z_PERSON'] [cat['HAS_S82XVI']==1]
z_best[cat['HAS_S82XVI']==1] = out
out2 = cat['TDVI_Z_PERSON']  [cat['HAS_TDVI']==1] 
z_best[cat['HAS_TDVI']==1]  = out2

cat_categories = cat['CLASS']
#out = cat['S82XVI_CLASS_PERSON'][cat['HAS_S82XVI']==1]
#cat_categories[cat['HAS_S82XVI']==1] = out
out = cat['TDVI_CLASS_PERSON'][cat['HAS_TDVI']==1]
cat_categories[cat['HAS_TDVI']==1] = out.strip()

print(survey)
for zmin, zmax in zip(zmins, zmaxs):
	print('=================================================')
	print(zmin, '<z<', zmax)
	selection = (z_best > zmin) & (z_best < zmax) & (cat['HAS_GOOD_Z']==1)
	all_categories = n.unique(cat_categories[selection])
	print(all_categories)
	for category in all_categories:
		ok_i = (cat_categories==category) & (selection)
		N_obj = len(ok_i.nonzero()[0])
		print(category, N_obj)
		if N_obj>0:
			name = survey+"_"+category+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
			print(name)
			DATA = n.transpose([ cat['PLATE'][ok_i], cat['MJD'][ok_i], cat['FIBERID'][ok_i], z_best[ok_i] ])
			n.savetxt(os.path.join(path_2_stack_lists, name), DATA)


# XXL lists
survey = 'XXL'
cat = cat_XXL

z_best = cat['Z']
out = cat['M16VI_Z_PERSON'] [cat['HAS_M16VI']==1]
z_best[cat['HAS_M16VI']==1] = out
out2 = cat['TDVI_Z_PERSON']  [cat['HAS_TDVI']==1] 
z_best[cat['HAS_TDVI']==1]  = out2

cat_categories = cat['CLASS']
out = cat['M16VI_CLASS_PERSON'][cat['HAS_M16VI']==1]
cat_categories[cat['HAS_M16VI']==1] = out
out = cat['TDVI_CLASS_PERSON'][cat['HAS_TDVI']==1]
cat_categories[cat['HAS_TDVI']==1] = out.strip()

print(survey)
for zmin, zmax in zip(zmins, zmaxs):
	print('=================================================')
	print(zmin, '<z<', zmax)
	selection = (z_best > zmin) & (z_best < zmax) & (cat['HAS_GOOD_Z']==1)
	all_categories = n.unique(cat_categories[selection])
	print(all_categories)
	for category in all_categories:
		ok_i = (cat_categories==category) & (selection)
		N_obj = len(ok_i.nonzero()[0])
		print(category, N_obj)
		if N_obj>0:
			name = survey+"_"+category+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
			print(name)
			DATA = n.transpose([ cat['PLATE'][ok_i], cat['MJD'][ok_i], cat['FIBERID'][ok_i], z_best[ok_i] ])
			n.savetxt(os.path.join(path_2_stack_lists, name), DATA)

#sys.exit()


# XMMSL2 lists
survey = 'XMMSL2'
cat = cat_XMMSL2
z_best = cat['Z_BEST']
z_conf = cat['CONF_BEST']
print(survey)
for zmin, zmax in zip(zmins, zmaxs):
	print('=================================================')
	print(zmin, '<z<', zmax)
	selection = (z_best > zmin) & (z_best < zmax) & (z_conf>=3)
	all_categories = n.unique(cat['CLASS_BEST'][selection])              
	for category in all_categories:
		ok_i = (cat['CLASS_BEST']==category) & (selection)
		N_obj = len(ok_i.nonzero()[0])
		print(category, N_obj)
		if N_obj>0:
			name = survey+"_"+category+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
			print(name)
			DATA = n.transpose([ cat['PLATE_BEST'][ok_i], cat['MJD_BEST'][ok_i], cat['FIBERID_BEST'][ok_i], z_best[ok_i] ])
			n.savetxt(os.path.join(path_2_stack_lists, name), DATA)


# 2RXS lists
survey = '2RXS'
cat = cat_2RXS
z_best = cat['Z_BEST']
z_conf = cat['CONF_BEST']
print(survey)
for zmin, zmax in zip(zmins, zmaxs):
	print('=================================================')
	print(zmin, '<z<', zmax)
	selection = (z_best > zmin) & (z_best < zmax) & (z_conf>=3)
	all_categories = n.unique(cat['CLASS_BEST'][selection])              
	for category in all_categories:
		ok_i = (cat['CLASS_BEST']==category) & (selection)
		N_obj = len(ok_i.nonzero()[0])
		print(category, N_obj)
		if N_obj>0:
			name = survey+"_"+category+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
			print(name)
			DATA = n.transpose([ cat['PLATE_BEST'][ok_i], cat['MJD_BEST'][ok_i], cat['FIBERID_BEST'][ok_i], z_best[ok_i] ])
			n.savetxt(os.path.join(path_2_stack_lists, name), DATA)

import glob

for zmin, zmax in zip(zmins, zmaxs):
	print('=================================================')
	name = "*_*_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
	all_ascii_lists = n.array(glob.glob(os.path.join(path_2_stack_lists, name)))
	all_ascii_lists.sort()
	if len(all_ascii_lists)>0:
		all_ascii_lists_survey = n.array([os.path.basename(el).split('_')[0] for el in all_ascii_lists])
		all_ascii_lists_category = n.array([os.path.basename(el).split('_')[1] for el in all_ascii_lists])
		all_categories = n.unique(all_ascii_lists_category)
		print(all_categories)
		for category in all_categories:
			ok_i = (all_ascii_lists_category==category) 
			files_2_join = " ".join(all_ascii_lists[ok_i])
			path_2_full_file = os.path.join(path_2_stack_lists, "full_"+category+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.asc' )
			command = "cat "+files_2_join+" > "+path_2_full_file
			print(command)
			os.system(command)

