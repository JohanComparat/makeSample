"""
python3.4 create_stack_lists_X_AGN.py
"""
from catalog_lib import *


data, high_conf, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, bl, nl, class_best = get_arrays_2rxs(path_2_cat_2RXS = path_2_cat_2RXS, EXI_ML_min = EXI_ML_min)

data_RXS = Catalog()
data_RXS.data = data
data_RXS.high_conf = high_conf
data_RXS.targeted = targeted
data_RXS.observed = observed
data_RXS.goodZ = goodZ
data_RXS.idZ = idZ
data_RXS.agnZ = agnZ
data_RXS.clusters = clusters
data_RXS.blazars_noZ = blazars_noZ
data_RXS.stars = stars
data_RXS.bl = bl
data_RXS.nl = nl
data_RXS.class_best = class_best

zmins = np.hstack((-0.01, 0.005, np.arange(0.1,3.0,0.1) ))
zmaxs = np.hstack((0.005, np.arange(0.,3.0,0.1)+0.5))
path_2_stack_lists = os.path.join(os.environ['HOME'], 'SDSS/stacks/X_AGN')
survey='ROSAT'
selection = data_RXS.clusters
category = 'clusterGAL' 

def write_stack_list(survey, selection, zmins, zmaxs, category = 'clusterGAL' ):
	plate = data_RXS.data['PLATE_BEST']          [selection]
	mjd = data_RXS.data['MJD_BEST']              [selection]
	fiberid = data_RXS.data['FIBERID_BEST']      [selection]
	class_best = data_RXS.data['merged_class']   [selection]
	subclass_star = data_RXS.data['XID_SUBCLASS'][selection]
	xid_star = data_RXS.data['XID_XID']          [selection]
	z_best = data_RXS.data['Z_BEST']             [selection]
	print(survey)
	for zmin, zmax in zip(zmins, zmaxs):
		print('=================================================')
		selection = (z_best > zmin) & (z_best < zmax) 
		N_obj = len(selection.nonzero()[0])
		print(zmin, '<z<', zmax, 'N', N_obj)
		if N_obj>10:
			name = survey+"_"+category+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
			print(name)
			DATA = np.transpose([ plate[selection], mjd[selection], fiberid[selection], z_best[selection] ])
			np.savetxt(os.path.join(path_2_stack_lists, name), DATA)

selection = data_RXS.clusters
UT = np.unique(data_RXS.data['merged_class'][selection])

category = 'CLUAGN'
selection = (data_RXS.clusters)&( (data_RXS.data['merged_class']=='BLAGN') | (data_RXS.data['merged_class']=='QSO') | (data_RXS.data['merged_class']=='NLAGN'))
max_z = np.max(data_RXS.data['Z_BEST']             [selection])
zmins = np.array([0.005])
zmaxs = np.array([0.5])
write_stack_list(survey, selection, zmins, zmaxs, category )

category = 'CLUGAL'
selection = (data_RXS.clusters)&(data_RXS.data['merged_class']=='GALAXY')
max_z = np.max(data_RXS.data['Z_BEST']             [selection])
zmins = np.array([0.005])
zmaxs = np.array([0.5])
write_stack_list(survey, selection, zmins, zmaxs, category )

sys.exit()

selection = data_RXS.clusters
category = 'clusterGAL' 
max_z = np.max(data_RXS.data['Z_BEST']             [selection])
zmins = np.hstack((-0.01, 0.005, np.arange(0.1,max_z,0.1) ))
zmaxs = np.hstack((0.005, np.arange(0.,max_z,0.1)+0.5))
write_stack_list(survey, selection, zmins, zmaxs, category )
zmaxs = np.hstack((0.005, np.arange(0.,max_z,0.1)+0.2))
write_stack_list(survey, selection, zmins, zmaxs, category )


selection = data_RXS.bl
category = 'AGNT1' 
max_z = np.max(data_RXS.data['Z_BEST']             [selection])
zmins = np.hstack((-0.01, 0.005, np.arange(0.1,max_z,0.1) ))
zmaxs = np.hstack((0.005, np.arange(0.,max_z,0.1)+0.5))
write_stack_list(survey, selection, zmins, zmaxs, category )
zmaxs = np.hstack((0.005, np.arange(0.,max_z,0.1)+0.2))
write_stack_list(survey, selection, zmins, zmaxs, category )


selection = data_RXS.nl
category = 'AGNT2' 
max_z = np.max(data_RXS.data['Z_BEST']             [selection])
zmins = np.hstack((-0.01, 0.005, np.arange(0.1,max_z,0.1) ))
zmaxs = np.hstack((0.005, np.arange(0.,max_z,0.1)+0.5))
write_stack_list(survey, selection, zmins, zmaxs, category )
zmaxs = np.hstack((0.005, np.arange(0.,max_z,0.1)+0.2))
write_stack_list(survey, selection, zmins, zmaxs, category )

sys.exit()

selection = data_RXS.stars
category = 'STARS' 
plate = data_RXS.data['PLATE_BEST']          [selection]
mjd = data_RXS.data['MJD_BEST']              [selection]
fiberid = data_RXS.data['FIBERID_BEST']      [selection]
class_best = data_RXS.data['merged_class']   [selection]
subclass_star = data_RXS.data['XID_SUBCLASS'][selection]
xid_star = data_RXS.data['XID_XID']          [selection]
z_best = data_RXS.data['Z_BEST']             [selection]

sub_CL = np.unique(subclass_star[(xid_star==1)])
print(survey)
for cl in sub_CL:
	print('=================================================')
	s1 = (xid_star==1)&(subclass_star==cl) 
	N_obj = len(s1.nonzero()[0])
	print(cl, 'N', N_obj)
	if cl=='CV/AM':
		cl_str='CVAM'
	elif cl=='CV/DN':
		cl_str='CVDN'
	else:
		cl_str=cl
	if N_obj>=10:
		name = survey+"_"+category+"_class_"+cl_str+'.ascii'
		print(name)
		DATA = np.transpose([ plate[s1], mjd[s1], fiberid[s1], z_best[s1] ])
		np.savetxt(os.path.join(path_2_stack_lists, name), DATA)





plate = data_RXS.data['PLATE_BEST']          [selection]
mjd = data_RXS.data['MJD_BEST']              [selection]
fiberid = data_RXS.data['FIBERID_BEST']      [selection]
class_best = data_RXS.data['merged_class']   [selection]
subclass_star = data_RXS.data['XID_SUBCLASS'][selection]
xid_star = data_RXS.data['XID_XID']          [selection]
z_best = data_RXS.data['Z_BEST']             [selection]

print(survey)
for zmin, zmax in zip(zmins, zmaxs):
	print('=================================================')
	selection = (z_best > zmin) & (z_best < zmax) 
	N_obj = len(selection.nonzero()[0])
	print(zmin, '<z<', zmax, 'N', N_obj)
	if N_obj>100:
		name = survey+"_"+category+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
		print(name)
		DATA = np.transpose([ plate[selection], mjd[selection], fiberid[selection], z_best[selection] ])
		np.savetxt(os.path.join(path_2_stack_lists, name), DATA)



survey='XMMSL2'

data, high_conf, targeted, observed, goodZ, idZ, agnZ, clusters, blazars_noZ, stars, bl, nl, class_best = get_arrays_xmmsl2(path_2_cat_xmmsl = path_2_cat_xmmsl, EXI_ML_min = EXI_ML_min)

data_XMM = Catalog()
data_XMM.data = data
data_XMM.high_conf = high_conf
data_XMM.targeted = targeted
data_XMM.observed = observed
data_XMM.goodZ = goodZ
data_XMM.idZ = idZ
data_XMM.agnZ = agnZ
data_XMM.clusters = clusters
data_XMM.blazars_noZ = blazars_noZ
data_XMM.stars = stars
data_XMM.bl = bl
data_XMM.nl = nl
data_XMM.class_best = class_best

sys.exit()
#zmins = np.array([-0.005])
#zmaxs = np.array([5.])

# 4 catalogs 
#path_2_XMMSL2 = os.path.join(path_2_cats, 'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR14_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')
#path_2_2RXS   = os.path.join(path_2_cats, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR14_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')
#path_2_XXL    = os.path.join(path_2_cats, 'Menzel16_spAll_XXL_plate.fits')
#path_2_S82X   = os.path.join(path_2_cats, 'LaMassa19_spAll_S82X_plate.fits')

#cat_XMMSL2 = fits.open(path_2_XMMSL2)[1].data
#cat_2RXS   = fits.open(path_2_2RXS  )[1].data
#cat_XXL    = fits.open(path_2_XXL   )[1].data
#cat_S82X   = fits.open(path_2_S82X  )[1].data

# 

# S82X lists
#survey = 'S82X'
#cat = cat_S82X

#z_best = cat['Z']
#out = cat['S82XVI_Z_PERSON'] [cat['HAS_S82XVI']==1]
#z_best[cat['HAS_S82XVI']==1] = out
#out2 = cat['TDVI_Z_PERSON']  [cat['HAS_TDVI']==1] 
#z_best[cat['HAS_TDVI']==1]  = out2

#cat_categories = cat['CLASS']
#out = cat['S82XVI_CLASS_PERSON'][cat['HAS_S82XVI']==1]
#cat_categories[cat['HAS_S82XVI']==1] = out
#out = cat['TDVI_CLASS_PERSON'][cat['HAS_TDVI']==1]
#cat_categories[cat['HAS_TDVI']==1] = out.strip()

# on laptop
#path_2_cats = os.path.join(os.environ['HOME'], 'data', 'spiders', 'agn')
# on servers
#path_2_cats = os.path.join(os.environ['HOME'], 'data1/SDSS/dr16/spiders')

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
	all_categories = np.unique(cat_categories[selection])
	print(all_categories)
	for category in all_categories:
		ok_i = (cat_categories==category) & (selection)
		N_obj = len(ok_i.nonzero()[0])
		print(category, N_obj)
		if N_obj>0:
			name = survey+"_"+category+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
			print(name)
			DATA = np.transpose([ cat['PLATE'][ok_i], cat['MJD'][ok_i], cat['FIBERID'][ok_i], z_best[ok_i] ])
			np.savetxt(os.path.join(path_2_stack_lists, name), DATA)

#sys.exit()


# XMMSL2 lists
survey = 'XMMSL2'
cat = cat_XMMSL2

z_best = cat['Z_BEST']
z_conf = cat['CONF_BEST']
cat_categories = cat['CLASS_BEST']

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

for zmin, zmax in zip(zmins[:3], zmaxs[:3]):
	print('=================================================')
	name = "*_*_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
	all_ascii_lists = n.array(glob.glob(os.path.join(path_2_stack_lists+'_v0', name)))
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


for zmin, zmax in zip(zmins, zmaxs):
	print('=================================================')
	name = "*_BLAGN_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
	all_ascii_lists1 = n.array(glob.glob(os.path.join(path_2_stack_lists+'_v0', name)))
	name = "*_QSO_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
	all_ascii_lists2 = n.array(glob.glob(os.path.join(path_2_stack_lists+'_v0', name)))
	all_ascii_lists = n.hstack((all_ascii_lists1, all_ascii_lists2))
	all_ascii_lists.sort()
	if len(all_ascii_lists)>0:
		files_2_join = " ".join(all_ascii_lists)
		path_2_full_file = os.path.join( path_2_stack_lists, "full_quasar_zmin_" + str(int(10*zmin)).zfill(2) + "_zmax_" + str(int(10*zmax)).zfill(2) + '.asc' )
		command = "cat " + files_2_join + " > " + path_2_full_file
		print(command)
		os.system(command)

