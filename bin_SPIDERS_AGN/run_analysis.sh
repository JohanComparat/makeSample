#!/bin/bash
cd ~/software/linux/makeSample/bin_SPIDERS_AGN

# specific eBOSS/SPIDERS AGN targets specZ catalog creation for clustering tests:
sh compile_specz_catalog_RASS_9028.sh
sh compile_specz_catalog_XMMSL_819.sh

# basic statistics and sky maps.

# section 2.
# section 2.1. 
# figure 1 :
python basic_stats_targets.py

# section 2.2. mask. inside / outside : 21279 / 9
# Table 1 and figure 2 :
python basic_stats_2RXS.py
# figures of completeness and success rate vs. flux and sky position
# python dr16-completeness-rass9028-xmmsl819.py 4 6p5
# python dr16-completeness-VAC-XMMSL2.py 4 6p5
# python dr16-completeness-VAC-2RXS.py 4 6p5
#
# python dr16-completeness-rass9028-xmmsl819.py 4 10
# python dr16-completeness-VAC-XMMSL2.py 4 10
# python dr16-completeness-VAC-2RXS.py 4 10

# section 3.
# table 1, 2:
# split by CONF_BEST
python dr16-completeness-VAC-2RXS.py 4 10
python dr16-completeness-VAC-2RXS.py 4 6p5

python dr16-completeness-VAC-XMMSL2.py 4 

# Figure 3
python visual-inspection-stats.py 6p5
python visual-inspection-stats.py 10
# # Figure 4, classification per type
# # python AGN-types.py 6p5

# figure 5.
# python plot_NZ-VAC.py 0.1 22.5 6.5
# python plot_NZ-VAC.py 0.2 22.5 6.5 
# python plot_NZ-VAC.py 0.0 22.5 6.5
# python plot_NZ-VAC.py 0.05 22.5 6.5
python plot_NZ-VAC.py 0.1 22.5 6.5
# python plot_NZ-VAC.py 0.15 22.5 6.5
# python plot_NZ-VAC.py 0.2 22.5 6.5

# 
# python plot_NZ-VAC.py 0.3 22.5 6.5
# python plot_NZ-VAC.py 0.4 22.5 6.5
# 
# python plot_NZ-VAC.py 0.45 22.5 6.5
# python plot_NZ-VAC.py 0.45 22.0 6.5
# python plot_NZ-VAC.py 0.45 21.5 6.5
# 
# python plot_NZ-VAC.py 0.5 22.5 6.5
# python plot_NZ-VAC.py 0.5 22.0 6.5
# python plot_NZ-VAC.py 0.5 21.5 6.5
# python plot_NZ-VAC.py 0.5 21.0 6.5
#                              6.5
# python plot_NZ-VAC.py 0.1 21.5 6.5
# python plot_NZ-VAC.py 0.2 21.5 6.5 
# python plot_NZ-VAC.py 0.3 21.5 6.5
# python plot_NZ-VAC.py 0.4 21.5 6.5
# python plot_NZ-VAC.py 0.5 21.5 6.5
#
# python plot_NZ-VAC.py 0.1 22.5 10.
# python plot_NZ-VAC.py 0.2 22.5 10.
# python plot_NZ-VAC.py 0.3 22.5 10.
# python plot_NZ-VAC.py 0.4 22.5 10.
# python plot_NZ-VAC.py 0.5 22.5 10.
# #                              10.
# python plot_NZ-VAC.py 0.1 21.5 10.
# python plot_NZ-VAC.py 0.2 21.5 10.
# python plot_NZ-VAC.py 0.3 21.5 10.
# python plot_NZ-VAC.py 0.4 21.5 10.
# python plot_NZ-VAC.py 0.5 21.5 10.
# 

# global X AGN spectroscopic catalogue construction
sh compile_specz_catalog_2RXS.sh
sh compile_specz_catalog_XMMSL2.sh

sh compile_random_catalog.sh

python create_gaia_mask_2RXS.py

# adds a boolean column if yes or no it is too close to a star
python apply_gaia_mask.py

python write_angular_clustering_catalogs.py
python compute_angular_clustering.py

python unwise-depth.py

# compute angular clustering up to a degree ?

# catalogs :
# /data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/2RXS/2RXS_AllWISE_catalog_paper_2017May26.fits.gz
# /data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/XMM/XMMSL2_AllWISE_catalog_paper_2017JUN09.fits.gz
# /data44s/eroAGN_WG_DATA/DATA/randoms/randoms.fits

# sub catalogs
# - full cat
# - cut abs(gal_lat)>20 degrees
# - cut on X-ray depth map
# - cut on wise depth map

# write ascii catalogs + param files here :
# /data36s/comparat/AGN_clustering/angular_clustering/

/data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/unwise/release/depth-mask/mask-0000m107.cat.fits
depth.list

UNWISE_DEPTH_CAT=/data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/unwise/release/depth-cat.fits

stilts tcat \
ifmt=fits in=@/data44s/eroAGN_WG_DATA/DATA/photometry/catalogs/unwise/release/depth.list \
omode=out out=$UNWISE_DEPTH_CAT \


# SDSS version = DR14
cd ~/software/linux/makeSample/bin_SPIDERS_AGN
# compares the AGN host mass function with the prediction from the mock catalog
run plot_Host_AGN_SMF_with_Bo16_model.py 0.001 0.4 25.

# creates the catalog to estimate clustering
python write_clustering_catalogs.py  0. 0.36
python write_clustering_catalogs.py  0.36 0.77
python write_clustering_catalogs.py  0.77 1.34

cd $OBS_REPO/SDSS/dr14/spiders/clustering_catalogs

cat clustering_agn_N_RL_0.005_N.data clustering_agn_N_RL_0.005_S.data > clustering_agn_RL_0.005.data  
cat clustering_agn_N_RL_0.005_N.random clustering_agn_N_RL_0.005_S.random > clustering_agn_RL_0.005.random

cat clustering_agn_N_RL_0.36_0.77_0.005_N.data  clustering_agn_N_RL_0.36_0.77_0.005_S.data > clustering_agn_RL_0.36_0.77_0.005.data

cat clustering_agn_N_RL_0.36_0.77_0.005_N.random clustering_agn_N_RL_0.36_0.77_0.005_S.random > clustering_agn_RL_0.36_0.77_0.005.random

cat clustering_agn_N_RL_0.77_1.34_0.005_N.data clustering_agn_N_RL_0.77_1.34_0.005_S.data > clustering_agn_RL_0.77_1.34_0.005.data

cat clustering_agn_N_RL_0.77_1.34_0.005_N.random clustering_agn_N_RL_0.77_1.34_0.005_S.random > clustering_agn_RL_0.77_1.34_0.005.random


cd $OBS_REPO/SDSS/dr14/spiders/clustering_measurements/

$DARKSIM_DIR/software/CUTE/CUTE/CUTE param_0_005_036_077.ini
$DARKSIM_DIR/software/CUTE/CUTE/CUTE param_0_005_077_134.ini
$DARKSIM_DIR/software/CUTE/CUTE/CUTE param_0_005.ini

cd ~/software/linux/makeSample/bin_SPIDERS_AGN
python plot_clustering.py