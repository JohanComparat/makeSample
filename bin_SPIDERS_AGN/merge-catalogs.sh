# Match with latest SDSS pipeline results, find 7208 counterparts

stilts tmatch2 \
in1=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26.fits.gz ifmt1=fits \
in2=/data37s/SDSS/catalogs/spAll_short_ZW04-v5_10_10.fits ifmt2=fits \
matcher=sky params="2" join=all1 find=best \
values1="ALLW_RA ALLW_DEC" values2="PLUG_RA PLUG_DEC" \
ocmd='addcol hp3 "healpixNestIndex( 3, ALLW_RA, ALLW_DEC )"' \
ocmd='addcol hp4 "healpixNestIndex( 4, ALLW_RA, ALLW_DEC )"' \
ocmd='addcol hp5 "healpixNestIndex( 5, ALLW_RA, ALLW_DEC )"' \
ocmd='addcol hp6 "healpixNestIndex( 6, ALLW_RA, ALLW_DEC )"' \
ocmd='addcol hp12 "healpixNestIndex( 12, ALLW_RA, ALLW_DEC )"' \
ocmd='addcol in_v5_10_10 "Separation>=0"' \
ocmd='delcols "Separation"' \
out=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_v5_10_10.fits 

stilts tmatch2 \
in1=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_v5_10_10.fits ifmt1=fits \
in2=/data37s/SDSS/catalogs/specObj-SDSS-dr14_short_ZW04.fits ifmt2=fits \
matcher=sky params="2" join=all1 find=best \
values1="ALLW_RA ALLW_DEC" values2="PLUG_RA PLUG_DEC" \
ocmd='addcol in_26 "Separation>=0"' \
ocmd='delcols "Separation"' \
out=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_v5_10_10_26.fits 


# then

# Match in 2 arcseconds swith tilts tool to the Veron & Veron catalog VII/258/vv10 to find 13,553 counterparts (10% of the targets)
# http://www.mpe.mpg.de/XraySurveys/2RXS_XMMSL2/ADD-ONs/2RXS_AllWISE_catalog_paper_2017May26.fits.gz
# http://cdsarc.u-strasbg.fr/viz-bin/Cat?VII/258

stilts tmatch2 \
in1=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_v5_10_10_26.fits ifmt1=fits \
in2=/data36s/comparat/veron-veron-13th-ed.fits ifmt2=fits \
icmd2='delcols "_RAJ2000 _DEJ2000 Name n_RAJ2000 FC"' \
matcher=sky params="2" join=all1 find=best \
values1="ALLW_RA ALLW_DEC" values2="RAJ2000 DEJ2000" \
ocmd='addcol in_veron "Separation>=0"' \
ocmd='delcols "Separation"' \
out=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_v5_10_10_26_VERON.fits 

# 7208 match v5_10_10
# 21721 match 26
# 13553 match Veron
# all matches: 21368


# create randoms in a healpix grid downsampling at N_with_z / N_targets. DEC>-20 is sufficient

# full sky randoms :

stilts tpipe \
in=/data37s/SDSS/catalogs/random-ra-dec.txt ifmt=ascii \
cmd='addcol hp3 "healpixNestIndex( 3, RA, DEC )"' \
cmd='addcol hp4 "healpixNestIndex( 4, RA, DEC )"' \
cmd='addcol hp5 "healpixNestIndex( 5, RA, DEC )"' \
cmd='addcol hp6 "healpixNestIndex( 6, RA, DEC )"' \
cmd='addcol hp12 "healpixNestIndex( 12, RA, DEC )"' \
omode=out out=/data37s/SDSS/catalogs/randoms.fits 

#### 2RXS MASK_2RXS_RATELIM
# /data36s/comparat/masks/2RXS_HPX_nside4096_all_ratelimit.fits.gz


stilts tmatch2 \
in1=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_v5_10_10_26_VERON.fits ifmt1=fits \
in2=/data36s/comparat/masks/2RXS_HPX_nside4096_all_ratelimit.fits ifmt2=fits \
matcher=exact join=all1 find=best \
values1="hp12" values2='$0' \
out=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_v5_10_10_26_VERON_2RXS_mask.fits 


stilts tmatch2 \
in1=/data37s/SDSS/catalogs/randoms.fits ifmt1=fits \
in2=/data36s/comparat/masks/2RXS_HPX_nside4096_all_ratelimit.fits ifmt2=fits \
matcher=exact join=all1 find=best \
values1="hp12" values2='$0' \
out=/data37s/SDSS/catalogs/randoms_2RXS_mask.fits 

# 
# compute tsr x ssr in hp12 and hp6 vs magnitude (wise) and X flux 
# ALLW_w1mpro
# ALLW_w1sigmpro
# 2RXS_SRC_FLUX
# 2RXS_SRC_FLUX
# set(hp6)
# make histogram per pixel of wise mag
# 
# 
# 
# 

import numpy as n
import astropy.io.fits as fits


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

import os

out_dir  = os.path.join(os.environ['HOME'], 'wwwDir/eRoMok/clustering/data/')

d=fits.open('/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_v5_10_10_26_VERON_2RXS_mask.fits')[1].data

spec = (d['in_veron']) | (d['in_26']) | (d['in_v5_10_10'])

hp3_list = n.array(list(set(d['hp3'][spec])))
hp3_list.sort()

hp4_list = n.array(list(set(d['hp4'][spec])))
hp4_list.sort()

# 4615  hp5 pixels with spectro
# 10607 hp6 pixels with spectro
w1_bins = n.hstack(( n.arange(0, 10., 5), n.arange(10, 21., 2.)))
x_w1 = (w1_bins[1:]+w1_bins[:-1])*0.5

def get_hist(pix, hp_name = 'hp3'):
	s=(d[hp_name]==pix)
	nall  = n.histogram(d['ALLW_w1mpro'][s], bins = w1_bins)[0]
	nspec = n.histogram(d['ALLW_w1mpro'][s & spec], bins = w1_bins)[0]
	return nspec, nall

hh = n.array([get_hist(pix) for pix in hp3_list ])




p.figure(1, (5,5))
p.axes([0.17,0.17,0.78,0.78])
for el in hh[:,1]:
	p.plot(x_w1, el, color='k', lw=0.1, rasterized=True)

for el in hh[:,0]:
	p.plot(x_w1, el, color='b', lw=0.1, rasterized=True)

p.xlabel('W1 mag')
p.ylabel('N/pixel ')
p.yscale('log')
#p.ylim((1e-4, 1))
p.xlim((4,20))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( out_dir,"w1_histogram.png"))
p.clf()


p.figure(1, (5,5))
p.axes([0.17,0.17,0.78,0.78])

ratio = []
for e1, e2 in zip(hh[:,0],hh[:,1]):
	den = e2
	den[e2==0]=1
	ratio.append(e1/den) 
	p.plot(x_w1, e1/den, color='k', lw=0.05, rasterized=True)

# deep fields :
sel = (hh[:,0].T[-2]>1)

ratio = n.array(ratio)
ratio_M = n.mean(ratio[sel], axis=0)
ratio_MD = n.median(ratio[sel], axis=0)
ratio_STD = n.std(ratio[sel], axis=0)

# p.errorbar(x_w1, ratio_M, yerr=ratio_STD, label='deep mean')
p.errorbar(x_w1, ratio_MD, yerr=ratio_STD, label='deep median')

# shallow fields :
sel = (hh[:,0].T[-2]<=1)

ratio = n.array(ratio)
ratio_M = n.mean(ratio[sel], axis=0)
ratio_MD = n.median(ratio[sel], axis=0)
ratio_STD = n.std(ratio[sel], axis=0)

# p.errorbar(x_w1, ratio_M, yerr=ratio_STD, label='shallow mean')
p.errorbar(x_w1, ratio_MD, yerr=ratio_STD, label='shallow median')

p.xlabel('W1 mag')
p.ylabel('N spec/N target ')
#p.yscale('log')
#p.ylim((1e-4, 1))
p.xlim((4,20))
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( out_dir,"w1_histogram_ratio.png"))
p.clf()



p.figure(1, (5,5))
p.axes([0.17,0.17,0.78,0.78])

# deep fields :
sel = (hh[:,0].T[-2]>1)
for pix in hp3_list[sel]:
	ss=(d['hp3']==pix)
	p.plot(d['ALLW_RA'][ss],d['ALLW_DEC'][ss],'k,', rasterized=True)

# shallow fields :
sel = (hh[:,0].T[-2]<=1)
for pix in hp3_list[sel]:
	ss=(d['hp3']==pix)
	p.plot(d['ALLW_RA'][ss],d['ALLW_DEC'][ss],'r,', rasterized=True)


p.xlabel('r.a.')
p.ylabel('dec')
p.grid()
p.legend(frameon=False, loc=0)
#p.title('data dr14')
p.savefig(os.path.join( out_dir,"w1_depth_spec_ra_dec.png"))
p.clf()









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

# match to /data37s/SDSS/catalogs/spAll-v5_10_10.fits

stilts tpipe \
in=/data37s/SDSS/26/catalogs/specObj-SDSS-dr14.fits \
cmd='delcols "SURVEY                     INSTRUMENT                  CHUNK                       PROGRAMNAME                 PLATERUN                    PLATEQUALITY             PLATESN2                   DEREDSN2                    LAMBDA_EFF                  BLUEFIBER                  ZOFFSET                     SNTURNOFF                  NTURNOFF                    SPECPRIMARY                 SPECSDSS                    SPECLEGACY                 SPECSEGUE                   SPECSEGUE1               SPECSEGUE2                 SPECBOSS                   BOSS_SPECOBJ_ID             SPECOBJID                   FLUXOBJID                   BESTOBJID                  TARGETOBJID                 PLATEID                     NSPECOBS                    FIRSTRELEASE               RUN2D                       RUN1D                      DESIGNID                    CX                         CY                          CZ                          XFOCAL                      YFOCAL                     SOURCETYPE                  TARGETTYPE                 THING_ID_TARGETING          THING_ID                   PRIMTARGET                  SECTARGET                   LEGACY_TARGET1              LEGACY_TARGET2             SPECIAL_TARGET1             SPECIAL_TARGET2            SEGUE1_TARGET1              SEGUE1_TARGET2             SEGUE2_TARGET1              SEGUE2_TARGET2              MARVELS_TARGET1             MARVELS_TARGET2            BOSS_TARGET1                BOSS_TARGET2               EBOSS_TARGET0               EBOSS_TARGET1              EBOSS_TARGET2               EBOSS_TARGET_ID             ANCILLARY_TARGET1           ANCILLARY_TARGET2          SPECTROGRAPHID             CHI68P                     FRACNSIGMA                  FRACNSIGHI                  FRACNSIGLO                 SPECTROFLUX                 SPECTROFLUX_IVAR           SPECTROSYNFLUX              SPECTROSYNFLUX_IVAR        SPECTROSKYFLUX              ANYANDMASK                 ANYORMASK                   SPEC1_G                    SPEC1_R                     SPEC1_I                    SPEC2_G                    SPEC2_R                    SPEC2_I                     ELODIE_FILENAME            ELODIE_OBJECT               ELODIE_SPTYPE              ELODIE_BV                   ELODIE_TEFF                ELODIE_LOGG                 ELODIE_FEH                 ELODIE_Z                    ELODIE_Z_ERR               ELODIE_Z_MODELERR           ELODIE_RCHI2               ELODIE_DOF                  Z_NOQSO                     Z_ERR_NOQSO                 ZWARNING_NOQSO             CLASS_NOQSO                 SUBCLASS_NOQSO             RCHI2DIFF_NOQSO             Z_PERSON                   CLASS_PERSON                Z_CONF_PERSON               COMMENTS_PERSON             CALIBFLUX                  CALIBFLUX_IVAR      "' \
cmd='select "(ZWARNING==0 || ZWARNING==4) && Z>Z_ERR && Z_ERR>0"' \
omode=out out=/data37s/SDSS/catalogs/specObj-SDSS-dr14_short_ZW04.fits 

stilts tmatch2 \
in1=/data37s/SDSS/catalogs/2RXS_AllWISE_catalog_paper_2017May26.fits ifmt1=fits \
in2=/data37s/SDSS/catalogs/spAll_short_ZW04-v5_10_10.fits ifmt2=fits \
matcher=sky params="1" join=all1 find=best suffix2=_DR14P \
values1="ALLW_RA ALLW_DEC" values2="PLUG_RA PLUG_DEC" \
ocmd='addcol hp5 "healpixNestIndex( 5, ALLW_RA, ALLW_DEC )"' \
ocmd='addcol hp6 "healpixNestIndex( 6, ALLW_RA, ALLW_DEC )"' \
ocmd='addcol in_v5_10_10 "Separation>=0"' \
ocmd='delcols "Separation"' \
out=/data37s/SDSS/catalogs/2RXS_AllWISE_catalog_paper_2017May26_DR14P.fits 

# tpipe 
# /data37s/SDSS/26/catalogs/specObj-SDSS-dr14_SNR.fits

#keepcols : 
     1 PLATE                      J
      2 TILE                       J
      3 MJD                        J
      4 FIBERID                    J
      5 RUN2D                      8A
      6 RUN1D                      8A
      7 OBJTYPE                    16A
      8 PLUG_RA                    D
      9 PLUG_DEC                   D
     10 CLASS                      6A
     11 SUBCLASS                   21A
     12 Z                          E
     13 Z_ERR                      E
     14 RCHI2                      E
     15 DOF                        J
     16 RCHI2DIFF                  E
     17 TFILE                      24A
     18 TCOLUMN                    10J
     19 NPOLY                      J
     20 THETA                      10E
     21 VDISP                      E
     22 VDISP_ERR                  E
     23 VDISPZ                     E
     24 VDISPZ_ERR                 E
     25 VDISPCHI2                  E
     26 VDISPNPIX                  E
     27 VDISPDOF                   J
     28 WAVEMIN                    E
     29 WAVEMAX                    E
     30 WCOVERAGE                  E
     31 ZWARNING                   J
     32 SN_MEDIAN                  5E
     33 SN_MEDIAN_ALL              E
     34 CHI68P                     E
     35 RA                         D
     36 DEC                        D
     37 CX                         D
     38 CY                         D
     39 CZ                         D
     40 RAERR                      D
     41 DECERR                     D
     42 L                          D
     43 B                          D
     44 OFFSETRA                   5E
     45 OFFSETDEC                  5E
     46 PSF_FWHM                   5E
     47 AIRMASS                    5E


     Column Name                Format     Dims       Units     TLMIN  TLMAX
      1 SURVEY                     6A
      2 INSTRUMENT                 4A
      3 CHUNK                      16A
      4 PROGRAMNAME                23A
      5 PLATERUN                   16A
      6 PLATEQUALITY               8A
      7 PLATESN2                   E
      8 DEREDSN2                   E
      9 LAMBDA_EFF                 E
     10 BLUEFIBER                  J
     11 ZOFFSET                    E
     12 SNTURNOFF                  E
     13 NTURNOFF                   J
     14 SPECPRIMARY                B
     15 SPECSDSS                   B
     16 SPECLEGACY                 B
     17 SPECSEGUE                  B
     18 SPECSEGUE1                 B
     19 SPECSEGUE2                 B
     20 SPECBOSS                   B
     21 BOSS_SPECOBJ_ID            J
     22 SPECOBJID                  22A
     23 FLUXOBJID                  19A
     24 BESTOBJID                  19A
     25 TARGETOBJID                22A
     26 PLATEID                    19A
     27 NSPECOBS                   I
     28 FIRSTRELEASE               3A
     29 RUN2D                      3A
     30 RUN1D                      A
     31 DESIGNID                   J
     32 CX                         D
     33 CY                         D
     34 CZ                         D
     35 XFOCAL                     E
     36 YFOCAL                     E
     37 SOURCETYPE                 19A
     38 TARGETTYPE                 8A
     39 THING_ID_TARGETING         K
     40 THING_ID                   J
     41 PRIMTARGET                 J
     42 SECTARGET                  J
     43 LEGACY_TARGET1             J
     44 LEGACY_TARGET2             J
     45 SPECIAL_TARGET1            K
     46 SPECIAL_TARGET2            K
     47 SEGUE1_TARGET1             J
     48 SEGUE1_TARGET2             J
     49 SEGUE2_TARGET1             J
     50 SEGUE2_TARGET2             J
     51 MARVELS_TARGET1            J
     52 MARVELS_TARGET2            J
     53 BOSS_TARGET1               K
     54 BOSS_TARGET2               K
     55 EBOSS_TARGET0              K
     56 EBOSS_TARGET1              K
     57 EBOSS_TARGET2              K
     58 EBOSS_TARGET_ID            K
     59 ANCILLARY_TARGET1          K
     60 ANCILLARY_TARGET2          K
     61 SPECTROGRAPHID             I
     62 PLATE                      J
     63 TILE                       J
     64 MJD                        J
     65 FIBERID                    J
     66 OBJID                      5J
     67 PLUG_RA                    D
     68 PLUG_DEC                   D
     69 CLASS                      6A
     70 SUBCLASS                   21A
     71 Z                          E
     72 Z_ERR                      E
     73 RCHI2                      E
     74 DOF                        J
     75 RCHI2DIFF                  E
     76 TFILE                      24A
     77 TCOLUMN                    10J
     78 NPOLY                      J
     79 THETA                      10E
     80 VDISP                      E
     81 VDISP_ERR                  E
     82 VDISPZ                     E
     83 VDISPZ_ERR                 E
     84 VDISPCHI2                  E
     85 VDISPNPIX                  E
     86 VDISPDOF                   J
     87 WAVEMIN                    E
     88 WAVEMAX                    E
     89 WCOVERAGE                  E
     90 ZWARNING                   J
     91 SN_MEDIAN_ALL              E
     92 SN_MEDIAN                  5E
     93 CHI68P                     E
     94 FRACNSIGMA                 10E
     95 FRACNSIGHI                 10E
     96 FRACNSIGLO                 10E
     97 SPECTROFLUX                5E
     98 SPECTROFLUX_IVAR           5E
     99 SPECTROSYNFLUX             5E
    100 SPECTROSYNFLUX_IVAR        5E
    101 SPECTROSKYFLUX             5E
    102 ANYANDMASK                 J
    103 ANYORMASK                  J
    104 SPEC1_G                    E
    105 SPEC1_R                    E
    106 SPEC1_I                    E
    107 SPEC2_G                    E
    108 SPEC2_R                    E
    109 SPEC2_I                    E
    110 ELODIE_FILENAME            25A
    111 ELODIE_OBJECT              21A
    112 ELODIE_SPTYPE              10A
    113 ELODIE_BV                  E
    114 ELODIE_TEFF                E
    115 ELODIE_LOGG                E
    116 ELODIE_FEH                 E
    117 ELODIE_Z                   E
    118 ELODIE_Z_ERR               E
    119 ELODIE_Z_MODELERR          E
    120 ELODIE_RCHI2               E
    121 ELODIE_DOF                 J
    122 Z_NOQSO                    E
    123 Z_ERR_NOQSO                E
    124 ZWARNING_NOQSO             J
    125 CLASS_NOQSO                A
    126 SUBCLASS_NOQSO             A
    127 RCHI2DIFF_NOQSO            E
    128 Z_PERSON                   E
    129 CLASS_PERSON               J
    130 Z_CONF_PERSON              J
    131 COMMENTS_PERSON            A
    132 CALIBFLUX                  5E
    133 CALIBFLUX_IVAR             5E
    134 PL                         D
    135 MJ                         D
    136 FI                         D
    137 SNR_ALL                    D
    138 SNR_32_35                  D
    139 SNR_35_39                  D
    140 SNR_39_41                  D
    141 SNR_41_55                  D
    142 SNR_55_68                  D
    143 SNR_68_74                  D
    144 SNR_74_93                  D
    145 Separation                 D
 
