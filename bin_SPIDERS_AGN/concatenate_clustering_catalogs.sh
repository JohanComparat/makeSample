#!/bin/bash/
cd $OBS_REPO/SDSS/dr14/spiders/clustering_catalogs
cat clustering_agn_N_RL_0.005_LX_440_N.data    clustering_agn_N_RL_0.005_LX_440_S.data    > clustering_agn_RL_0.005_LX_440.data       
cat clustering_agn_N_RL_0.005_LX_440_N.random  clustering_agn_N_RL_0.005_LX_440_S.random  > clustering_agn_RL_0.005_LX_440.random     
cat clustering_agn_N_RL_0.005_LX_443_N.data    clustering_agn_N_RL_0.005_LX_443_S.data    > clustering_agn_RL_0.005_LX_443.data   
cat clustering_agn_N_RL_0.005_LX_443_N.random  clustering_agn_N_RL_0.005_LX_443_S.random  > clustering_agn_RL_0.005_LX_443.random 
cat clustering_agn_N_RL_0.005_LX_446_N.data    clustering_agn_N_RL_0.005_LX_446_S.data    > clustering_agn_RL_0.005_LX_446.data   
cat clustering_agn_N_RL_0.005_LX_446_N.random  clustering_agn_N_RL_0.005_LX_446_S.random  > clustering_agn_RL_0.005_LX_446.random 

