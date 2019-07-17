#!/bin/bash
cd ~/software/linux/makeSample/bin_SPIDERS_AGN

sh compile_specz_catalog_RASS_9028.sh
sh compile_specz_catalog_XMMSL_819.sh

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