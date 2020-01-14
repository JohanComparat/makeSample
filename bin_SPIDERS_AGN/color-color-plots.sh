#!/bin/bash

topcat -stilts plot2plane \
   xpix=600 ypix=600 \
   xlabel='SDSS MODELMAG g-r' ylabel='SDSS MODELMAG r-i' grid=true fontsize=14 fontstyle=serif \
    fontweight=bold \
   xmin=-0.5 xmax=2 ymin=-0.5 ymax=2 \
   legend=true legpos=0.0,1.0 \
   in=/home/comparat/data/AGN_clustering/catalogs/2RXS/SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1.fits \
    x=SDSS_MODELMAG_g-SDSS_MODELMAG_r y=SDSS_MODELMAG_r-SDSS_MODELMAG_i shading=auto \
   layer_1=Mark \
      shape_1=filled_diamond size_1=2 color_1=grey \
      leglabel_1='All counterparts' \
   layer_2=Mark \
      icmd_2='select is_bl' \
      size_2=4 color_2=cyan \
      leglabel_2='type 1 AGN' \
   layer_3=Mark \
      icmd_3='select is_nl' \
      size_3=4 color_3=black \
      leglabel_3='type 2 AGN' \
   layer_4=Mark \
      icmd_4='select is_blazar' \
      shape_4=filled_square size_4=3 color_4=grey \
      leglabel_4='blazar' \
   layer_5=Mark \
      icmd_5='select is_star' \
      size_5=4 color_5=magenta \
      leglabel_5='star' \
   layer_6=Mark \
      icmd_6='select "is_cluster && (merged_class==\"QSO\" || merged_class==\"BLAGN\")"' \
      shape_6=filled_triangle_up size_6=3 color_6=orange \
      leglabel_6='QSO in cluster' \
   layer_7=Mark \
      icmd_7='select "is_cluster && merged_class==\"GALAXY\""' \
      shape_7=filled_triangle_down size_7=4 color_7=red \
      leglabel_7='galaxy in cluster' \
   omode=out out=$GIT_MAKESAMPLE/figures/agn/figures_VAC/color-color-gr-ri-stilts.png

   
topcat -stilts plot2plane \
   xpix=600 ypix=600 \
   xlabel=W2-W3 ylabel=W1-W2 grid=true fontsize=14 fontstyle=serif fontweight=bold \
   xmin=-0.5 xmax=5 ymin=-0.5 ymax=2 \
   legend=true legpos=0.0,1.0 \
   in=/home/comparat/data/AGN_clustering/catalogs/2RXS/SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1.fits \
    x=ALLW_W2mag-ALLW_W3mag y=ALLW_W1mag-ALLW_W2mag shading=auto \
   layer_1=Mark \
      shape_1=filled_diamond size_1=2 color_1=grey \
      leglabel_1='All counterparts' \
   layer_2=Mark \
      icmd_2='select is_bl' \
      size_2=4 color_2=cyan \
      leglabel_2='type 1 AGN' \
   layer_3=Mark \
      icmd_3='select is_nl' \
      size_3=4 color_3=black \
      leglabel_3='type 2 AGN' \
   layer_4=Mark \
      icmd_4='select is_blazar' \
      shape_4=filled_square size_4=3 color_4=grey \
      leglabel_4='blazar' \
   layer_5=Mark \
      icmd_5='select is_star' \
      size_5=4 color_5=magenta \
      leglabel_5='star' \
   layer_6=Mark \
      icmd_6='select "is_cluster && (merged_class==\"QSO\" || merged_class==\"BLAGN\")"' \
      shape_6=filled_triangle_up size_6=3 color_6=orange \
      leglabel_6='QSO in cluster' \
   layer_7=Mark \
      icmd_7='select "is_cluster && merged_class==\"GALAXY\""' \
      shape_7=filled_triangle_down size_7=4 color_7=red \
      leglabel_7='galaxy in cluster' \
   omode=out out=$GIT_MAKESAMPLE/figures/agn/figures_VAC/color-color-w1w2-w2w3-stilts.png

topcat -stilts plot2plane \
   xpix=600 ypix=600 \
   xlabel='\log_{10}(F_X [erg/cm^2/s])' ylabel=W1 grid=true texttype=latex fontsize=14 fontstyle=serif \
    fontweight=bold \
   xmin=-13 xmax=-10.5 ymin=4 ymax=18 \
   legend=true legpos=1.0,1.0 \
   in=/home/comparat/data/AGN_clustering/catalogs/2RXS/SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE_MaxBCG_REDMAPPER_SPIDERSCODEX_XID_classifications_LXtype1.fits \
    x='log10(RXS_SRC_FLUX)' y=ALLW_W1mag shading=auto \
   layer_1=Mark \
      shape_1=filled_diamond size_1=2 color_1=grey \
      leglabel_1='All\; counterparts' \
   layer_2=Mark \
      icmd_2='select is_bl' \
      size_2=4 color_2=cyan \
      leglabel_2='type\; 1\; AGN' \
   layer_3=Mark \
      icmd_3='select is_nl' \
      size_3=4 color_3=black \
      leglabel_3='type\; 2\; AGN' \
   layer_4=Mark \
      icmd_4='select is_blazar' \
      shape_4=filled_square size_4=3 color_4=grey \
      leglabel_4='blazar' \
   layer_5=Mark \
      icmd_5='select is_star' \
      size_5=4 color_5=magenta \
      leglabel_5='star' \
   layer_6=Mark \
      icmd_6='select "is_cluster && (merged_class==\"QSO\" || merged_class==\"BLAGN\")"' \
      shape_6=filled_triangle_up size_6=3 color_6=orange \
      leglabel_6='QSO\; in \; cluster' \
   layer_7=Mark \
      icmd_7='select "is_cluster && merged_class==\"GALAXY\""' \
      shape_7=filled_triangle_down size_7=4 color_7=red \
      leglabel_7='galaxy\; in\; cluster' \
   layer_8=Function \
      fexpr_8='-1.65*x-8.8' color_8=black thick_8=3 dash_8=3,3 \
      leglabel_8='-1.65\log_{10}(F_X)-8.8' \
   omode=out out=$GIT_MAKESAMPLE/figures/agn/figures_VAC/color-color-w1-FX-stilts.png

