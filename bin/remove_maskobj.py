#!/usr/bin/env python

import astropy.io.fits as fits
import subprocess
import os
import sys
from os.path import join
import pymangle as mangle
import numpy as np

infits  = sys.argv[1]
tagra   = sys.argv[2]
tagdec  = sys.argv[3]
outfits = sys.argv[4]
masklist= sys.argv[5]
depthMask = sys.argv[6]

def usage():
    print 'Usage:'
    print 'remove_maskobj.py infits tagra tagdec outfits masklist'
    print ''
    print ' selects in infile unmasked objects, and print this in outfile'
    print ' - infits   (string): input .fits file'
    print ' - tagra    (string): ra field name in infile'
    print ' - tagdec   (string): dec field name in infile'
    print ' - outfits  (string): output .fits file'
    print ' - masklist (string): mangle ply masks list, comma separated (e.g.,"rykoff,Tycho20Vmag10,Tycho210Vmag11,Tycho211Vmag115")'
    return

# masks
maskdir    = join(os.environ['LEGACYSURVEYTS_DIR'],'masks')
masklist   = masklist.split(',')
## to be completed if more masks desired
dict_mask  = {}
dict_mask['rykoff']          = join(maskdir, 'bright_object_mask_rykoff_pix.ply')
dict_mask['Tycho20Vmag10']   = join(maskdir, 'tycho2mask-0Vmag10.pol')
dict_mask['Tycho210Vmag11']  = join(maskdir, 'tycho2mask-10Vmag11.pol')
dict_mask['Tycho211Vmag115'] = join(maskdir, 'tycho2mask-11Vmag115.pol')


# adding mask info
hdu     = fits.open(infits)
cols    = hdu[1].columns
data    = hdu[1].data
ra      = data[tagra]
dec     = data[tagdec]

ismask  = np.zeros(len(ra))
for mask in masklist:
    # we read the masks
    mng         = mangle.Mangle(dict_mask[mask])
    # we retrieve the mask indexes
    polyid      = mng.polyid(ra,dec)
    # we mark masked objects
    tmp         = (polyid!=-1)
    ismask[tmp] = ismask[tmp] + 1

# depth mask in decals
if depthMask=='NGC':
	gL = 62.79716079 
	rL = 30.05661087
	zL = 11.

if depthMask=='SGC':
	gL = 62.79716079 
	rL = 30.05661087
	zL = 12.75  

if depthMask=='none':
	gL = 0. 
	rL = 0.
	zL = 0.

sel = (data['depth_ivar_g_wBug']<=gL) & (data['depth_ivar_r_wBug']<=rL) & (data['depth_ivar_z_wBug']<=zL) 
ismask[sel]=ismask[sel]+1

if infits.split('/')[-1]=="elg_190_ngc.fits" or infits.split('/')[-1]=="random-sweep1.fits":
	FP = (110<ra)&(ra<180)&(10<dec)&(dec<35)
	out = (FP==False)
	ismask[out] = ismask[out]+1

#NGC (115<ra<175, 15<dec<30, 830 deg2)
if infits.split('/')[-1]=="elg_240_sgc.fits" or infits.split('/')[-1]=="random-sweep0.fits" or infits.split('/')[-1]=="random-sweep3.fits" :
	FP = ((316.5<ra)&(ra<=360)&(-2<dec)&(dec<2)) | ((ra>=0)&(ra<45)&(-5<dec)&(dec<5))
	out = (FP==False)
	ismask[out] = ismask[out]+1

# unmasked obj.
tmp     = (ismask==0)

# writing outfits with unmasked obj. only
subprocess.call('cp '+infits+' '+outfits, shell=True)
fits.update(outfits,data[tmp],1)
 
