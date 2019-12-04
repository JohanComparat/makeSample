import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse

spec_dir = join(os.environ['HOME'],"SDSS/stacks/X_AGN")
stack_dir = join(os.environ['GIT_MAKESAMPLE'],"data/X_AGN_stacks")

file_list = n.array(glob.glob(os.path.join(stack_dir, 'ROSAT_CLU*.ascii')))

def stack_it(specList ):
	outfile = join(spec_dir, os.path.basename(specList)[:-5]+".stack")
	print(outfile)
	test_D = n.loadtxt(specList, unpack=True)
	print(len(test_D[0]))
	if len(test_D[0])>10:
		stack=sse.SpectraStackingEBOSS(specList, outfile)
		stack.createStackMatrix()
		stack.stackSpectra()

for file_input in file_list:
	stack_it(file_input)


