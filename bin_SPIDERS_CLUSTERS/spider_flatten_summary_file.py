import astropy.io.fits as fits
import os, sys
import time
import numpy as n
t0 = time.time()

#dir = os.path.join(os.environ['HOME'], 'hegcl', 'SPIDERS', 'codex_firefly_matching_2018_mar_22')

#hd = fits.open(os.path.join(dir, 'CODEX-DR14-MergedSpectroscopicCatalogue_2018-02-05.fits'))[1].data


dir2 = os.path.join(os.environ['HOME'], 'hegcl', 'SPIDERS')
dir2 = os.path.join(os.environ['HOME'], 'data', 'spiders', 'cluster')

#hd = fits.open(os.path.join(dir, 'validatedclusters_catalogue_2018-04-27_version_round123-v1_Xmass123-v1.fits'))[1].data
#hd = fits.open(os.path.join(dir, 'validatedclusters_catalogue_2018-12-04_version_round1-v1_Xmass1-v1.fits.gz'))[1].data
hd = fits.open(os.path.join(dir2, 'mastercatalogue_FINAL_CODEXID.fits'))[1].data

DATA = []
for el in hd :
  print( el['CLUS_ID'] )
  line_start = n.array([
    el['CLUS_ID'], 
    el['CODEX'], 
    el['KIND'], 
    el['NCOMPONENT'], 
    el['RA'], 
    el['DEC'], 
    el['RA_OPT'], 
    el['DEC_OPT'], 
    el['LAMBDA_CHISQ'], 
    el['LAMBDA_CHISQ_OPT'], 
    el['Z_LAMBDA'], 
    el['Z_LAMBDA_ERR'], 
    el['OBSSTATUS'], 
    el['NSUBMITTED'], 
    el['NTILED'], 
    el['NHASZ'], 
    el['NISNEWEBOSS'], 
    el['NISNEWSEQUELS'], 
    el['NISNEWSPIDERS'], 
    el['NISNEWSPIDERS_RASS_CLUS'], 
    el['NISNEWSPIDERS_XCLASS_CLUS'], 
    el['NISNEWSEQUELS_RASS_CLUS'], 
    el['NMEM'], 
    el['CLUZSPEC'], 
    el['CLUZSPECBOOT'], 
    el['CLUZSPECRUEL'], 
    el['CLUVDISP_GAP'], 
    el['CLUVDISPERR_GAP'], 
    el['CLUVDISP_BWT'], 
    el['CLUVDISPERR_BWT'], 
    el['CLUVDISPTYPE'], 
    el['CLUVDISPBEST'], 
    el['CLUVDISPBESTERR'], 
    el['ADDCOMMENT'], 
    el['DAZSPEC'], 
    el['NOKZ'], 
    el['NMEMBERS'], 
    el['SCREEN_CLU_RA_W'], 
    el['SCREEN_CLU_DEC_W'], 
    el['SCREEN_CLUZSPEC'], 
    el['SCREEN_CLUZSPEC_ERR'], 
    el['SCREEN_CLUZSPEC_SPREAD'], 
    el['SCREEN_CLUVDISP_GAP'], 
    el['SCREEN_CLUVDISP_BWT'], 
    el['SCREEN_CLUVDISP_BEST'], 
    el['SCREEN_NMEMBERS_W'], 
    el['STATUS'], 
    el['NINSPECTORS'], 
    el['NVALID'], 
    el['M200C'], 
    el['EM200C'], 
    el['LX0124'], 
    el['ELX'], 
    el['LC'], 
    el['R200C_DEG'], 
    el['KT'], 
    el['EKT'], 
    el['EZ'], 
    el['RAPP_R500'], 
    el['ACOR'], 
    el['KCOR'], 
    el['PSFCOR'], 
    el['FLUX052'], 
    el['EFLUX052'], 
    el['SIGMA'], 
    el['NH'] 
    ])
  for ii, plate in enumerate(el['ALLPLATE']):
    if plate != -99999 and plate > 0  :
      line_part2 = n.array([
    el[ 'ALLZ'][ii], 
    el[ 'ALLZ_ERR'][ii], 
    el[ 'ALLZWARNING'][ii], 
    el[ 'ALLZ_NOQSO'][ii], 
    el[ 'ALLZ_ERR_NOQSO'][ii], 
    el[ 'ALLZWARNING_NOQSO'][ii], 
    el[ 'ALLPRIORITIES_SEQUELS'][ii], 
    el[ 'ALLPRIORITIES_SPIDERS'][ii], 
    el[ 'ALLFIBER2MAG_I'][ii], 
    el[ 'ALLCMODELMAG_I'][ii], 
    el[ 'ALLPROBABILITIES'][ii], 
    el[ 'ALLRA'][ii], 
    el[ 'ALLDEC'][ii], 
    el[ 'ALLBCGFLAG'][ii], 
    el[ 'ALLPLATE'][ii], 
    el[ 'ALLMJD'][ii], 
    el[ 'ALLFIBERID'][ii], 
    el[ 'ISMEMBER'][ii], 
    el['SCREEN_ISMEMBER_W'][ii] 
    ])
      new_line = n.hstack((line_start, line_part2))
      DATA.append(new_line) 

header1 = " CLUS_ID    CODEX     KIND    COMPONENT     RA     DEC     RA_OPT     DEC_OPT    LAMBDA_CHISQ     LAMBDA_CHISQ_OPT     Z_LAMBDA     Z_LAMBDA_ERR     OBSSTATUS     NSUBMITTED     NTILED     NHASZ     NISNEWEBOSS     NISNEWSEQUELS    NISNEWSPIDERS     NISNEWSPIDERS_RASS_CLUS    NISNEWSPIDERS_XCLASS_CLUS     NISNEWSEQUELS_RASS_CLUS     NMEM     CLUZSPEC     CLUZSPECBOOT     CLUZSPECRUEL     CLUVDISP_GAP     CLUVDISPERR_GAP     CLUVDISP_BWT     CLUVDISPERR_BWT     CLUVDISPTYPE     CLUVDISPBEST     CLUVDISPBESTERR    ADDCOMMENT     DAZSPEC     NOKZ     NMEMBERS     SCREEN_CLU_RA_W     SCREEN_CLU_DEC_W     SCREEN_CLUZSPEC     SCREEN_CLUZSPEC_ERR     SCREEN_CLUZSPEC_SPREAD     SCREEN_CLUVDISP_GAP     SCREEN_CLUVDISP_BWT     SCREEN_CLUVDISP_BEST     SCREEN_NMEMBERS_W     STATUS     NINSPECTORS     NVALID     M200C     EM200C     LX0124     ELX     LC     R200C_DEG     KT    EKT     EZ     RAPP_R500     ACOR     KCOR     PSFCOR     FLUX052     EFLUX052     SIGMA    NH   IDLSPEC1D_Z   IDLSPEC1D_Z_ERR   IDLSPEC1D_ZWARNING   IDLSPEC1D_Z_NOQSO   IDLSPEC1D_Z_ERR_NOQSO   IDLSPEC1D_ZWARNING_NOQSO  PRIORITIES_SEQUELS   PRIORITIES_SPIDERS   FIBER2MAG_I   CMODELMAG_I  PROBABILITIES   RA_GAL   DEC_GAL   BCGFLAG   PLATE   MJD   FIBERID   ISMEMBER  SCREEN_ISMEMBER_W "

header = ", ".join(header1.split())

#n.savetxt(os.path.join(dir, 'CODEX-DR14-MergedSpectroscopicCatalogue_2018-02-05-flat.csv'), DATA, header = header, fmt='%s', delimiter=',')

out_file = os.path.join(dir2, 'mastercatalogue_FINAL_CODEXID-flat.csv')
n.savetxt(out_file, DATA, header = header, fmt='%s', delimiter=',')

sys.exit()

out_file_fits = os.path.join(dir2, 'mastercatalogue_FINAL_CODEXID-flat.fits')
cmd = """stilts tpipe in="""+out_file+""" ifmt=csv omode=out ofmt=fits out="""+out_file_fits
print(cmd)
os.system(cmd)

# then does the matching on the DS machines, spALL is 14G !!

spAll_file = '/home/comparat/data2/firefly/v1_1_0/v5_13_0/catalogs/sdss_eboss_firefly-dr16.fits'
out_file_fits_all = os.path.join(dir2, 'mastercatalogue_FINAL_CODEXID-flat_DR16_firefly2.fits')


c0 = """stilts tmatch2 ifmt1=fits ifmt2=fits in1="""+out_file_fits
c1 = """ in2="""+spAll_file+"""  omode=out out="""+out_file_fits_all
c2 = """ matcher=3d values1='PLATE MJD FIBERID' values2='PLATE MJD FIBERID' params=0.00001 join=all1  suffix1=_spiders suffix2=_DR16 """

cmd = c0+c1+c2
print(cmd)

#stilts tmatch2 ifmt1=fits ifmt2=fits in1=/home/comparat/hegcl/SPIDERS/mastercatalogue_FINAL_CODEXID-flat.fits in2=/home/comparat/data2/firefly/v1_1_0/v5_13_0/catalogs/sdss_eboss_firefly-dr16.fits  omode=out out=/home/comparat/hegcl/SPIDERS/mastercatalogue_FINAL_CODEXID-flat_DR16_firefly.fits matcher=3d values1='PLATE MJD FIBERID' values2='PLATE MJD FIBERID' params=0.00001 join=all1  suffix1=_spiders suffix2=_DR16 

