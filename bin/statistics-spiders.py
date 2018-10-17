import os, sys
import glob
import numpy as n
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as p
import astropy.io.fits as fits
plot_dir = os.path.join(os.environ['OBS_REPO'], 'SDSS', 'dr16')

## observed in the year 4 sample
#data = fits.open(os.path.join(plot_dir, 'spAll-v5_11_0.fits'))[1].data

#nl = lambda sel : len(sel.nonzero()[0])

#good_z = (  data['ZWARNING']==0 ) 

#eboss_qso_0 = ( (data['EBOSS_TARGET0'] & ( 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 ) ) > 0 )
#agn_0 = ( (data['EBOSS_TARGET0'] & ( 20 | 22 ) ) > 0 )
#clu_0 = ( (data['EBOSS_TARGET0'] & ( 21 | 23 ) ) > 0 )

#eboss_qso_1 = ( (data['EBOSS_TARGET1'] & ( 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 ) ) > 0 )
#s82x_1 = ( (data['EBOSS_TARGET1'] & ( 32 | 33 | 34 ) ) > 0 )

#agn_2 = ( (data['EBOSS_TARGET2'] & ( 0 | 2 | 4 ) ) > 0 )
#clu_2 = ( (data['EBOSS_TARGET2'] & ( 1 | 3 | 5 ) ) > 0 )
#tdss_qso = ( (data['EBOSS_TARGET2'] & ( 21 | 23 | 24 | 25 | 27 | 34 | 35 | 36 | 37 | 38 ) ) > 0 )
#s82x_2 = ( (data['EBOSS_TARGET2'] & ( 50 | 51 | 52 | 53 | 54 | 55 | 56 | 57 | 58 | 59 | 60 | 61 | 62 ) ) > 0 )

#x_agn = (agn_0 | agn_2)
#x_clu = (clu_0 | clu_2)
#eboss_qso = (eboss_qso_0 | eboss_qso_1)
#s82_x = (s82x_1 | s82x_2 )

#print('X AGN'    , nl(x_agn)     , nl(good_z & x_agn    ), int(100.*nl(good_z & x_agn)    *1./nl(x_agn)    )  )
#print('X CLUSTER', nl(x_clu)     , nl(good_z & x_clu    ), int(100.*nl(good_z & x_clu)    *1./nl(x_clu)    )  )
#print('S82X'     , nl(s82_x)     , nl(good_z & s82_x    ), int(100.*nl(good_z & s82_x)    *1./nl(s82_x)    )  )
#print('TDSS QSO' , nl(tdss_qso)  , nl(good_z & tdss_qso ), int(100.*nl(good_z & tdss_qso) *1./nl(tdss_qso) )  )
#print('eboss QSO', nl(eboss_qso) , nl(good_z & eboss_qso), int(100.*nl(good_z & eboss_qso)*1./nl(eboss_qso))  )

#print('X AGN & CLU'    , nl(x_agn & x_clu)   )

#redshift = data['Z']

#zbins = n.arange(0,3,0.1)
#p.figure(1, (8,5))
#p.hist(redshift[good_z & eboss_qso], bins = zbins, label='eboss_qso', histtype='step', lw=3)
#p.hist(redshift[good_z & s82_x    ], bins = zbins, label='s82_x    ', histtype='step', lw=2)
#p.hist(redshift[good_z & tdss_qso ], bins = zbins, label='tdss_qso ', histtype='step', lw=1, )
#p.hist(redshift[good_z & x_agn    ], bins = zbins, label='x_agn    ', histtype='step', lw=3)
#p.hist(redshift[good_z & x_clu    ], bins = zbins, label='x_clu    ', histtype='step', lw=2)
#p.xlabel('redshift')
#p.ylabel('number')
#p.yscale('log')
#p.grid()
#p.legend(frameon=False)
#p.savefig(os.path.join(plot_dir, 'NZ.png'))
#p.clf()

# futur tiles
tile_dir = os.path.join(os.environ['OBS_REPO'], 'SDSS', 'tilelist')
final_tiles = n.array(glob.glob(os.path.join(tile_dir, 'final-eboss?.fits')))
final_tiles.sort()

p.figure(1, (10,7))

for tile_f in final_tiles:
  hd = fits.open(tile_f)
  Npt = len(hd[1].data['DEC'])
  rd = n.random.random(Npt)
  sel = (rd < 2000./Npt)
  p.plot(hd[1].data['RA'][sel], hd[1].data['DEC'][sel], rasterized=True, marker='.', linestyle='None', label=os.path.basename(tile_f)[6:-5])
  print(tile_f, Npt, len(hd[1].data['RA'][sel]))

p.xlim((0,360))
p.ylim((-10,90))
p.xlabel('ra')
p.ylabel('dec')
p.grid()
p.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=4, fancybox=True, shadow=True)
p.savefig(os.path.join(plot_dir, 'chunks-eboss.png'))
p.clf()

# futur tiles
tile_dir = os.path.join(os.environ['OBS_REPO'], 'SDSS', 'tilelist')
final_tiles = n.array(glob.glob(os.path.join(tile_dir, 'final-eboss1?.fits')))
final_tiles.sort()

p.figure(1, (10,7))

for tile_f in final_tiles:
  hd = fits.open(tile_f)
  Npt = len(hd[1].data['DEC'])
  rd = n.random.random(Npt)
  sel = (rd < 2000./Npt)
  p.plot(hd[1].data['RA'][sel], hd[1].data['DEC'][sel], rasterized=True, marker='.', linestyle='None', label=os.path.basename(tile_f)[6:-5])
  print(tile_f, Npt, len(hd[1].data['RA'][sel]))

p.xlim((0,360))
p.ylim((-10,90))
p.xlabel('ra')
p.ylabel('dec')
p.grid()
p.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=4, fancybox=True, shadow=True)
p.savefig(os.path.join(plot_dir, 'chunks-eboss1.png'))
p.clf()

# futur tiles
tile_dir = os.path.join(os.environ['OBS_REPO'], 'SDSS', 'tilelist')
final_tiles = n.array(glob.glob(os.path.join(tile_dir, 'final-eboss2?.fits')))
final_tiles.sort()

p.figure(1, (10,7))

for tile_f in final_tiles:
  hd = fits.open(tile_f)
  Npt = len(hd[1].data['DEC'])
  rd = n.random.random(Npt)
  sel = (rd < 2000./Npt)
  p.plot(hd[1].data['RA'][sel], hd[1].data['DEC'][sel], rasterized=True, marker='.', linestyle='None', label=os.path.basename(tile_f)[6:-5])
  print(tile_f, Npt, len(hd[1].data['RA'][sel]))

p.xlim((0,360))
p.ylim((-10,90))
p.xlabel('ra')
p.ylabel('dec')
p.grid()
p.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=4, fancybox=True, shadow=True)
p.savefig(os.path.join(plot_dir, 'chunks-eboss2.png'))
p.clf()


# futur tiles
tile_dir = os.path.join(os.environ['OBS_REPO'], 'SDSS', 'tilelist')
final_tiles = n.array(glob.glob(os.path.join(tile_dir, 'final-eboss???.fits')))
final_tiles.sort()

p.figure(1, (10,7))

for tile_f in final_tiles:
  hd = fits.open(tile_f)
  Npt = len(hd[1].data['DEC'])
  rd = n.random.random(Npt)
  sel = (rd < 2000./Npt)
  p.plot(hd[1].data['RA'][sel], hd[1].data['DEC'][sel], rasterized=True, marker='.', linestyle='None', label=os.path.basename(tile_f)[6:-5])
  print(tile_f, Npt, len(hd[1].data['RA'][sel]))

p.xlim((0,360))
p.ylim((-10,90))
p.xlabel('ra')
p.ylabel('dec')
p.grid()
p.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=4, fancybox=True, shadow=True)
p.savefig(os.path.join(plot_dir, 'chunks-eboss3.png'))
p.clf()
