
# Match in 2 arcseconds with topcat cds tool to the Veron & Veron catalog VII/258/vv10 to find 13,403 counterparts
stilts tmatch2 \
in1=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26.fits.gz ifmt1=fits \
in2=/data36s/comparat/veron-veron-13th-ed.fits ifmt2=fits \
icmd2='delcols "_RAJ2000 _DEJ2000 Name n_RAJ2000 FC"' \
matcher=sky params="2" join=all1 find=best \
values1="ALLW_RA ALLW_DEC" values2="RAJ2000 DEJ2000" \
ocmd='addcol in_veron "Separation>=0"' \
out=/data36s/comparat/SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_VERON.fits 



INFO
====

for finding the counterparts to 
 - the 2RXS sources, we adopted an area of 30,434.6 sqdeg.
 - the XMMSL2 sources, we adopted an area of 25,565 sqdeg.

MISSING MASK : WISE MASK and ROSAT FIELD MASK !!
 
Data
====

full target catalog :
http://www.mpe.mpg.de/XraySurveys/2RXS_XMMSL2/ADD-ONs/2RXS_AllWISE_catalog_paper_2017May26.fits.gz
with 
132,254 targets
ALLW_RA, ALLW_DEC

Match in 2 arcseconds with topcat cds tool to the Veron & Veron catalog VII/258/vv10
to find 13,403 counterparts

Match to VAC_spiders_2RXS_DR14.fits
to find at separation 0 with ALLW_RA, ALLW_DEC: 
9,055 / 9,073 counterparts

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
