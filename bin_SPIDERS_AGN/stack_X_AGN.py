import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse

spec_dir = join(os.environ['HOME'],"SDSS/stacks/X_AGN")

file_list = n.array(glob.glob(os.path.join(spec_dir, 'full_*_zmin_00_zmax_50.asc')))

def stack_it(specList ):
 outfile = join(spec_dir, os.path.basename(specList)[:-4]+".stack")
 print(outfile)
 test_D = n.loadtxt(specList, unpack=True)
 print(len(test_D[0]))
 if len(test_D[0])>10:
  stack=sse.SpectraStackingEBOSS(specList, outfile, l_start=2.8, l_end=4.1 )
  stack.createStackMatrix()
  stack.stackSpectra()

for file_input in file_list:
 stack_it(file_input)


