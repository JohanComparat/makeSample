import os, sys

out_dir = '/data36s/comparat/AGN_clustering/angular_clustering/'

def compute_clustering(name, out_dir):
	out_name = os.path.join(out_dir , name )
	f=open(out_name+'.ini','w')
	f.write('data_filename= '+out_name+'.data \n')    
	f.write('random_filename= '+out_name+'.random \n')    
	f.write('input_format= 2 \n')
	f.write('output_filename= '+out_name+'.wtheta \n')
	f.write('corr_type= angular \n')
	f.write('omega_M= 0.307 \n')
	f.write('omega_L= 0.693 \n')
	f.write('w= -1 \n')
	f.write('log_bin= 10 \n')
	f.write('dim1_min_logbin= 0.001 \n')
	f.write('dim1_max= 1. \n')
	f.write('dim1_nbin= 60 \n')
	f.write('dim2_max= 160. \n')
	f.write('dim2_nbin= 40 \n')
	f.write('dim3_min= 0.00 \n')
	f.write('dim3_max= 3. \n')
	f.write('dim3_nbin= 1 \n')
	f.write('use_pm= 0 \n')
	f.write('n_pix_sph= 1000000 \n')
	f.close()
	os.system("~/darksim/software/CUTE/CUTE/CUTE "+out_name+'.ini')

compute_clustering('2RXS_AllWISE_catalog_paper_2017May26_X_GAL', out_dir)
compute_clustering('XMMSL2_AllWISE_catalog_paper_2017JUN09_X_GAL', out_dir)


def plot_results(name, out_dir):
	out_name = os.path.join(out_dir , name )

	fig = p.figure(2, (10,8))
	# wedge plot
	fig.add_subplot(2,2,1)

	ra, dec, z, w = np.loadtxt(out_name + '.data', unpack=True)
	p.plot(ra, dec, marker=',', color='b', alpha=0.1, rasterized = True, ls='None')
	p.xlabel('R.A. [deg]')
	p.ylabel('Dec. [deg]')
	p.title('data')

	fig.add_subplot(2,2,2)
	ra, dec, z, w = np.loadtxt(out_name + '.random', unpack=True)
	p.plot(ra, dec, marker=',', color='b', alpha=0.1, rasterized = True, ls='None')
	p.xlabel('R.A. [deg]')
	p.ylabel('Dec. [deg]')
	p.title('randoms')

	fig.add_subplot(2,2,3)
	# nz
	DATA = np.loadtxt(out_name + '.wtheta', unpack=True)
	p.errorbar(DATA[0], DATA[1], yerr=DATA[1]*DATA[2]**(-0.5), rasterized = True)
	p.xlabel('theta [deg]')
	p.ylabel('w(theta) ')
	p.xscale('log')
	p.yscale('log')
	p.grid()
	
	p.savefig(out_name+".png")
	p.clf()

plot_results('XMMSL2_AllWISE_catalog_paper_2017JUN09_X_GAL', out_dir)
plot_results('2RXS_AllWISE_catalog_paper_2017May26_X_GAL', out_dir)


