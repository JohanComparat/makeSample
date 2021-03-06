import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse

spec_dir = join(os.environ['HOME'],"SDSS/stacks/SPIDERS_C_GAL_z")

#file_list = n.array(glob.glob(os.path.join(spec_dir, 'spiders_C_GAL_??_z_??.ascii')))
file_list = n.array(glob.glob(os.path.join(spec_dir, '*.csv')))

file_list = n.array(glob.glob(os.path.join('/home/comparat/software/linux/makeSample', 'data/allobjects_mem075_alllx_package', '*.csv')))
file_list.sort()

def stack_it(specList ):
 outfile = join(spec_dir, os.path.basename(specList)[:-4]+".stack")
 print(outfile)
 #test_D = n.loadtxt(specList, unpack=True)
 #print(len(test_D[0]))
 #if len(test_D[0])>10:
 stack=sse.SpectraStackingEBOSS(specList, outfile, l_start=3., l_end=4.1, csv_input=True )
 stack.createStackMatrix()
 stack.stackSpectra()

#file_name = os.path.join(spec_dir, 'clusterCGAL_allz.ascii')
#stack_it(file_name)
#file_name = os.path.join(spec_dir, 'clusterCGAL_allz_inner.ascii')
#stack_it(file_name)

for file_input in file_list:
 stack_it(file_input)



