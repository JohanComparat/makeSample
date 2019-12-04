#!/bin/bash
cd ~/software/linux/makeSample/bin_SPIDERS_AGN

# specific eBOSS/SPIDERS AGN targets specZ catalog creation for clustering tests:
sh compile_specz_catalog_RASS_9028.sh
sh compile_specz_catalog_XMMSL_819.sh

# basic statistics and sky maps.

# section 2.
# section 2.1. 
# pseudo, figure 1 :
# python 000_basic_stats_targets.py

# section 2.2. mask. inside / outside : 21279 / 9
# Table 1 :
python 001_basic_stats_2RXS.py
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
python 002_dr16-completeness-VAC-2RXS.py 4 10
python 002_dr16-completeness-VAC-2RXS.py 4 6p5
python 002_dr16-completeness-VAC-XMMSL2.py 4 

python 002_dr16-completeness-FULLBOSS.py 4 10


# Figure 2, 3
python 003_visual-inspection-stats.py 6p5
python 003_visual-inspection-stats.py 10
# # Figure 4, classification per type
# # python AGN-types.py 6p5

python 004_NZ-VAC.py 0.1 22.5 6.5

# figure 5.
# python 004_NZ-VAC.py 0.1 22.5 6.5
# python 004_NZ-VAC.py 0.2 22.5 6.5 
# python 004_NZ-VAC.py 0.0 22.5 6.5
# python 004_NZ-VAC.py 0.05 22.5 6.5
# python 004_NZ-VAC.py 0.15 22.5 6.5
# python 004_NZ-VAC.py 0.2 22.5 6.5

# 
# python 004_NZ-VAC.py 0.3 22.5 6.5
# python 004_NZ-VAC.py 0.4 22.5 6.5
# 
# python 004_NZ-VAC.py 0.45 22.5 6.5
# python 004_NZ-VAC.py 0.45 22.0 6.5
# python 004_NZ-VAC.py 0.45 21.5 6.5
# 
# python 004_NZ-VAC.py 0.5 22.5 6.5
# python 004_NZ-VAC.py 0.5 22.0 6.5
# python 004_NZ-VAC.py 0.5 21.5 6.5
# python 004_NZ-VAC.py 0.5 21.0 6.5
#                              6.5
# python 004_NZ-VAC.py 0.1 21.5 6.5
# python 004_NZ-VAC.py 0.2 21.5 6.5 
# python 004_NZ-VAC.py 0.3 21.5 6.5
# python 004_NZ-VAC.py 0.4 21.5 6.5
# python 004_NZ-VAC.py 0.5 21.5 6.5
#
# python 004_NZ-VAC.py 0.1 22.5 10.
# python 004_NZ-VAC.py 0.2 22.5 10.
# python 004_NZ-VAC.py 0.3 22.5 10.
# python 004_NZ-VAC.py 0.4 22.5 10.
# python 004_NZ-VAC.py 0.5 22.5 10.
# #                              10.
# python 004_NZ-VAC.py 0.1 21.5 10.
# python 004_NZ-VAC.py 0.2 21.5 10.
# python 004_NZ-VAC.py 0.3 21.5 10.
# python 004_NZ-VAC.py 0.4 21.5 10.
# python 004_NZ-VAC.py 0.5 21.5 10.
# 

# figure LX vs. Z
python plot_LX_Z-VAC.py

topcat -stilts plot2plane \
   xpix=500 ypix=500 \
   xlog=true xlabel=redshift \
    ylabel='LX (soft), log10(Luminosity/[erg/s]) ' grid=true \
    texttype=latex fontsize=15 fontweight=bold \
   xmin=0.01 xmax=5.5 ymin=39 ymax=47.5 \
   legend=true legpos=0.0,1.0 \
   shading=auto \
   layer_1=Mark \
      in_1=/home/comparat/data/eRoMok/mocks/eRosita_eRASS8_with_photometry.fits \
      x_1=redshift_R y_1=LX_soft \
      color_1=grey \
      leglabel_1=eROSITA \
   layer_2=Mark \
      in_2=/home/comparat/software/makeSample/data/lxz/deep_fields_specz.dat \
       ifmt_2=ASCII \
      x_2=z y_2=logLX \
      shape_2=cross \
   layer_3=Mark \
      in_3=/home/comparat/software/makeSample/data/lxz/XMM_XXL_lz.dat \
       ifmt_3=ASCII \
      x_3=Z y_3=logLUM_SOFT \
      shape_3=cross \
      leglabel_3='deep\; surveys' \
   layer_4=Mark \
      in_4=/home/comparat/data/AGN_clustering/catalogs/2RXS/SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1.fits \
      x_4=Z_BEST y_4='log10(DC_2RXS_LX)' \
      size_4=2 color_4=blue \
      leglabel_4=2RXS \
   layer_5=Mark \
      in_5=/home/comparat/data/AGN_clustering/catalogs/XMMSL2/SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1_not2RXS.fits \
      x_5=Z_BEST y_5='log10(DC_XMMSL2_LX)' \
      size_5=2 color_5=yellow \
      leglabel_5=XMMSL2 \
   legseq=_1,_3,_4,_5 \
   omode=out out=$GIT_MAKESAMPLE/figures/agn/figures_VAC/LX_redshift.png

topcat -stilts plot2plane \
   xpix=500 ypix=300 \
   xlog=true ylog=true xlabel=redshift ylabel=Counts grid=true texttype=latex \
    fontsize=15 fontweight=bold \
   xmin=0.01 xmax=5.5 ymin=10 ymax=1000000 \
   legend=false \
   binsize=-12 barform=steps thick=4 \
   layer_1=Histogram \
      in_1=/home/comparat/data/eRoMok/mocks/eRosita_eRASS8_with_photometry.fits \
      x_1=redshift_R \
      color_1=grey \
      leglabel_1=eROSITA \
   layer_2=Histogram \
      in_2=/home/comparat/software/makeSample/data/lxz/XMM_XXL_lz_and_deep_fields.fits \
      x_2=z \
      leglabel_2='deep\; surveys' \
   layer_3=Histogram \
      in_3=/home/comparat/data/AGN_clustering/catalogs/2RXS/SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1.fits \
       icmd_3='select DC_2RXS_LX>0' \
      x_3=Z_BEST \
      color_3=blue \
      leglabel_3=2RXS \
   layer_4=Histogram \
      in_4=/home/comparat/data/AGN_clustering/catalogs/XMMSL2/SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1_not2RXS.fits \
       icmd_4='select DC_XMMSL2_LX>0' \
      x_4=Z_BEST \
      color_4=yellow \
      leglabel_4=XMMSL2 \
      omode=out out=$GIT_MAKESAMPLE/figures/agn/figures_VAC/LX_redshift_histZ.png

topcat -stilts plot2plane \
   xpix=500 ypix=300 \
   xlog=true ylog=true xlabel= ylabel=Counts \
    grid=true texttype=latex fontsize=15 fontweight=bold \
   xmin=39 xmax=47.5 ymin=10 ymax=1000000 \
   legend=false \
   binsize=-14 barform=steps thick=4 \
   layer_1=Histogram \
      in_1=/home/comparat/data/eRoMok/mocks/eRosita_eRASS8_with_photometry.fits \
      x_1=LX_soft \
      color_1=grey \
   layer_2=Histogram \
      in_2=/home/comparat/software/makeSample/data/lxz/XMM_XXL_lz_and_deep_fields.fits \
      x_2=logLX \
   layer_3=Histogram \
      in_3=/home/comparat/data/AGN_clustering/catalogs/2RXS/SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1.fits \
       icmd_3='select DC_2RXS_LX>0' \
      x_3='log10(DC_2RXS_LX)' \
      color_3=blue \
   layer_4=Histogram \
      in_4=/home/comparat/data/AGN_clustering/catalogs/XMMSL2/SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1_not2RXS.fits \
       icmd_4='select DC_XMMSL2_LX>0' \
      x_4='log10(DC_XMMSL2_LX)' \
      color_4=yellow \
      omode=out out=$GIT_MAKESAMPLE/figures/agn/figures_VAC/LX_redshift_histLX.png
      
python plot_color-color-mag-categories.py

EXIT

# stack spectrA according to J Wolf's classification
python stack_X_AGN.py
python plot_stacks.py
python plot_stacks_EV1.py

# https://github.com/JohanComparat/makeSample/blob/master/bin_SPIDERS_AGN/create_stack_lists_X_AGN.py
# https://github.com/JohanComparat/makeSample/blob/master/bin_SPIDERS_AGN/stack_X_AGN.py
# https://github.com/JohanComparat/makeSample/blob/master/bin_SPIDERS_AGN/plot_stacks.py
# https://github.com/JohanComparat/pySU/blob/master/galaxy/python/SpectraStackingEBOSS.py


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