
# Match in 2 arcseconds swith tilts tool to the Veron & Veron catalog VII/258/vv10 to find 13,553 counterparts (10% of the targets)
# http://www.mpe.mpg.de/XraySurveys/2RXS_XMMSL2/ADD-ONs/2RXS_AllWISE_catalog_paper_2017May26.fits.gz
# http://cdsarc.u-strasbg.fr/viz-bin/Cat?VII/258

stilts tmatch2 \
in1=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26.fits.gz ifmt1=fits \
in2=/data36s/comparat/veron-veron-13th-ed.fits ifmt2=fits \
icmd2='delcols "_RAJ2000 _DEJ2000 Name n_RAJ2000 FC"' \
matcher=sky params="2" join=all1 find=best \
values1="ALLW_RA ALLW_DEC" values2="RAJ2000 DEJ2000" \
ocmd='addcol in_veron "Separation>=0"' \
ocmd='delcols "Separation"' \
out=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_VERON.fits 

# Match in 2 arcseconds with stilts tool to the VAC_spiders_2RXS_DR14 to find 9,055 counterparts (7% of the targets)
# https://data.sdss.org/sas/dr14/eboss/spiders/analysis/VAC_spiders_2RXS_DR14.fits

stilts tmatch2 \
in1=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_VERON.fits ifmt1=fits \
in2=/data36s/comparat/SDSS/dr14/spiders/analysis/VAC_spiders_2RXS_DR14.fits ifmt2=fits \
icmd2='delcols "RXS_IAU_NAME  RXS_ExiML               RXS_Cts                 RXS_e_Cts               RXS_CRate               RXS_e_CRate             RXS_ExpTime             RXS_RAJ2000             RXS_DEJ2000             RXS_Ext                 RXS_e_Ext               RXS_ExtML               RXS_RADEC_ERR           RXS_NEAR_BRIGHT_STAR    RXS_NEAR_BRIGHT_GAL     LOGGALNH                RXS_SRC_FLUX            RXS_SRC_FLUX_ERR        NWAY_bias               NWAY_p_any              NWAY_p_i                ALLW_designation        ALLW_sigra              ALLW_sigdec             ALLW_RADECERR           ALLW_w1mpro             ALLW_w1sigmpro          ALLW_w1snr              ALLW_w2mpro             ALLW_w2sigmpro          ALLW_w2snr              ALLW_w3mpro             ALLW_w3sigmpro          ALLW_w3snr              ALLW_w4mpro             ALLW_w4sigmpro          ALLW_w4snr              ALLW_cc_flags           ALLW_ext_flg            ALLW_var_flg            ALLW_r_2mass            ALLW_n_2mass            ALLW_j_m_2mass          ALLW_j_msig_2mass       ALLW_h_m_2mass          ALLW_h_msig_2mass       ALLW_k_m_2mass          ALLW_k_msig_2mass       NWAY_dist_XRAY_ALLW     GAIA_DR1_ra               GAIA_DR1_dec              GAIA_DR1_source_id        GAIA_DR1_ref_epoch        GAIA_DR1_ra_error         GAIA_DR1_dec_error        GAIA_DR1_parallax         GAIA_DR1_parallax_error   GAIA_DR1_pmra             GAIA_DR1_pmra_error       GAIA_DR1_pmdec            GAIA_DR1_pmdec_error      GAIA_DR1_phot_g_mean_flux GAIA_ALLW_angDist GAIA_DR1_phot_g_mean_flux_error"' \
matcher=sky params="2" join=all1 find=best suffix1="" suffix2=_DR14 \
values1="ALLW_RA ALLW_DEC" values2="ALLW_RA ALLW_DEC" \
ocmd='addcol in_DR14 "Separation>=0"' \
ocmd='delcols "ALLW_RA_DR14 ALLW_DEC_DR14 Separation"' \
out=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_VERON_DR14.fits 

# Match in 2 arcseconds with stilts tool to the VAC_spiders_XMMSL_DR14 to find  counterparts (% of the targets)
# https://data.sdss.org/sas/dr14/eboss/spiders/analysis/VAC_spiders_XMMSL_DR14.fits       

stilts tmatch2 \
in1=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_VERON_DR14.fits ifmt1=fits \
in2=/data36s/comparat/SDSS/dr14/spiders/analysis/VAC_spiders_XMMSL_DR14.fits ifmt2=fits \
icmd2='delcols "GAIA_DR1_ra               GAIA_DR1_dec              GAIA_DR1_source_id        GAIA_DR1_ref_epoch        GAIA_DR1_ra_error         GAIA_DR1_dec_error        GAIA_DR1_parallax         GAIA_DR1_parallax_error   GAIA_DR1_pmra             GAIA_DR1_pmra_error       GAIA_DR1_pmdec            GAIA_DR1_pmdec_error      GAIA_DR1_phot_g_mean_flux GAIA_ALLW_angDist GAIA_DR1_phot_g_mean_flux_error"' \
matcher=sky params="2" join=all1 find=best suffix1="" suffix2=_XMMSL \
values1="ALLW_RA ALLW_DEC" values2="ALLW_RA ALLW_DEC" \
ocmd='addcol in_XMMSL "Separation>=0"' \
ocmd='delcols "ALLW_RA_DR14 ALLW_DEC_DR14 Separation"' \
out=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_VERON_DR14_XMMSL.fits 

INFO
====

for finding the counterparts to 
 - the 2RXS sources, we adopted an area of 30,434.6 sqdeg.
 - the XMMSL2 sources, we adopted an area of 25,565 sqdeg.

MISSING MASK : WISE MASK and ROSAT FIELD MASK !!
 
Data
====

Match with specObj-dr14.fits
within 1 arcsecond to PLUG_RA, PLUG_DEC
to find 17,178 counterparts

add column
hpID_4096
healpixNestIndex(12, ALLW_RA, ALLW_DEC)

match to depth maps from 2RXS : 2RXS_HPX_nside4096_all_ratelimit.fits.gz
match in integer (exact)
hpID_4096 and $0

Save here the catalog :
/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask.fits
Spectroscopic results are available for 20,152 targets.

Mask for DR14 eBOSS footprint 
cd /data36s/comparat/SDSS/LSS/QSO/v1.9f4$ 
ls 
mask-QSO-N-eboss_v1.9f4.ply  mask-QSO-S-eboss_v1.9f4.ply

are copied here : /data36s/comparat/SDSS/dr14/spiders/masks/

cd /data36s/comparat/SDSS/dr14/spiders
python apply_dr14_ngc_sgc_mask.py 

creates a NGC and a SGC catalog
/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask_DR14areaNGC.fits
/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask_DR14areaSGC.fits

Randoms
=======

Use the randoms created for the eBOSS QSO
/data36s/comparat/SDSS/LSS/QSO/v1.9f4/
eboss_v1.9f4-QSO-N-eboss_v1.9f4.ran.fits       eboss_v1.9f4-QSO-S-eboss_v1.9f4.ran.fits      

add hpID_4096 :
healpixNestIndex(12, ALLW_RA, ALLW_DEC)

match to depth maps from 2RXS : 2RXS_HPX_nside4096_all_ratelimit.fits.gz
match in integer (exact)
hpID_4096 and $0

saved here :
/data36s/comparat/SDSS/dr14/spiders/randoms/

Finding the right depth cut
===========================
And writes the catalogs in ascii

python ks-test-depth-value.py

first of all angular KS test 

KS-test in ra and dec between 
 - data 0<z<0.2 and flux > 2e-13
 - randoms

for increasing values of 'MASK_2RXS_RATELIM' acceptance that mask both data and randoms.
