#!/bin/bash

topcat -stilts plot2sky \
   xpix=913 ypix=590 \
   projection=car \
   legend=true legpos=1.0,0.0 \
   in=/home/comparat/data/AGN_clustering/XMMSL2_AllWISE_catalog_paper_2017JUN09_v5_13_0_sdss_26_VERON_2QZ_MASKED.fits \
    lon=ALLW_RA lat=ALLW_DEC shading=auto size=2 \
   layer_1=Mark \
      color_1=grey \
      leglabel_1='3: All' \
   layer_2=Mark \
      icmd_2='select in_veron' \
      color_2=yellow \
      leglabel_2='3: in_veron' \
   layer_3=Mark \
      icmd_3='select in_SDSS_26' \
      color_3=blue \
      leglabel_3='3: in_SDSS_26' \
   layer_4=Mark \
      icmd_4='select in_BOSS_v5_13_0' \
      leglabel_4='3: in_BOSS_v5_13_0' \
   layer_5=Mark \
      icmd_5='select in_2QZ' \
      color_5=magenta \
      leglabel_5='3: in_2QZ' 

topcat -stilts plot2sky \
   xpix=913 ypix=590 \
   projection=car \
   legend=true legpos=1.0,0.0 \
   in=/home/comparat/data/AGN_clustering/2RXS_AllWISE_catalog_paper_2017May26_v5_13_0_sdss_26_VERON_2QZ_MASKED.fits \
    lon=ALLW_RA lat=ALLW_DEC shading=auto size=2 \
   layer_1=Mark \
      color_1=light_grey \
      leglabel_1='1: All' \
   layer_2=Mark \
      icmd_2='select in_veron' \
      color_2=yellow \
      leglabel_2='1: in_veron' \
   layer_3=Mark \
      icmd_3='select in_SDSS_26' \
      color_3=blue \
      leglabel_3='1: in_SDSS_26' \
   layer_4=Mark \
      icmd_4='select in_BOSS_v5_13_0' \
      leglabel_4='1: in_BOSS_v5_13_0' \
   layer_5=Mark \
      icmd_5='select in_2QZ' \
      color_5=magenta \
      leglabel_5='1: in_2QZ' 

