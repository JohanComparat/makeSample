#!/bin/bash/

python remove_maskobj.py \
$DECALS_DATA_DIR/catalogs/elg_240_sgc.fits \
ra \
dec \
$DECALS_DATA_DIR/catalogs/elg_240_sgc.masked.fits \
rykoff,Tycho20Vmag10,Tycho210Vmag11,Tycho211Vmag115 \
SGC \

python remove_maskobj.py \
$DECALS_DATA_DIR/catalogs/ngc_dr3bis/elg190ngcDR3B.fits \
ra \
dec \
$DECALS_DATA_DIR/catalogs/ngc_dr3bis/elg190ngcDR3B.masked.fits \
rykoff,Tycho20Vmag10,Tycho210Vmag11,Tycho211Vmag115 \
NGC \

python remove_maskobj.py \
$DECALS_DATA_DIR/catalogs/ngc_dr3/elg190ngcDR3.fits \
ra \
dec \
$DECALS_DATA_DIR/catalogs/ngc_dr3/elg190ngcDR3.masked.fits \
rykoff,Tycho20Vmag10,Tycho210Vmag11,Tycho211Vmag115 \
NGC \

python remove_maskobj.py \
$DECALS_RANDOMS_DIR/random-ngc.fits \
ra \
dec \
$DECALS_RANDOMS_DIR/random-ngc.masked.fits \
rykoff,Tycho20Vmag10,Tycho210Vmag11,Tycho211Vmag115 \
NGC \

python remove_maskobj.py \
$DECALS_DATA_DIR/randoms/random-sweep0.fits \
ra \
dec \
$DECALS_DATA_DIR/randoms/random-sweep0.maskedLow.fits \
rykoff,Tycho20Vmag10,Tycho210Vmag11,Tycho211Vmag115 \
SGC \

python remove_maskobj.py \
$DECALS_DATA_DIR/randoms/random-sweep3.fits \
ra \
dec \
$DECALS_DATA_DIR/randoms/random-sweep3.maskedLow.fits \
rykoff,Tycho20Vmag10,Tycho210Vmag11,Tycho211Vmag115 \
SGC \
