"""
Postprocessing of firefly's outputs.

1. measures chi2 and ndof for each fit

2. creates a catalog per list of spFly files, adds chi2 and ndof.

"""
import time
t0t=time.time()
from os.path import join
import sys, os
import numpy as n
import glob 
import astropy.io.fits as fits
from scipy.interpolate import interp1d
from astropy.table import Table, Column

# order of the hdus in the model files
#hdu_header_prefix = n.array(["Chabrier_MILES_", "Salpeter_UVMILES_" , "Kroupa_MILES_"   , "Chabrier_ELODIE_", "Salpeter_ELODIE_", "Kroupa_ELODIE_"  , "Chabrier_STELIB_", "Salpeter_STELIB_", "Kroupa_STELIB_" ])

path_2_spec_dir = join(os.environ['HOME'],"SDSS/stacks/SPIDERS_C_GAL")
path_2_spec_dir = os.path.join(os.environ['HOME'], 'software/makeSample/data/allobjects_mem075_alllx_package')

path_2_fly_dir = join(path_2_spec_dir, "firefly")

file_list = n.array(glob.glob(os.path.join(path_2_fly_dir, 'spFly-*.stack')))
baseNames = n.array([ os.path.basename(bn)[6:] for bn in file_list ])

hdu_header_prefix = n.array(["Chabrier_MILES_", "Chabrier_ELODIE_" ])

dV=-9999
z_name = 'Z_NOQSO'

path_2_out_file = join( path_2_fly_dir, "spFlyAll.fits")

# parse baseNames

title_2 = n.array([ baseName[12:-6] for baseName in baseNames ])
title_2_sp = n.array([ title.split('_') for title in title_2 ])

#title = str( n.round( float(title_2_sp[0]),3 ) ) + r'$<\log_{10}(\theta/\theta_{200c})<$' + str( n.round( float(title_2_sp[2] ),3 ) )

#sys.exit()

redshift_dict_min = {'z0-0.15' : 0.00, 'z0.15-0.3' :0.15, 'z0.3-0.45':0.30, 'z0.45+': 0.45 } 
redshift_dict_max = {'z0-0.15' : 0.15, 'z0.15-0.3' :0.30, 'z0.3-0.45':0.45, 'z0.45+': 0.70 } 
bin_dict_max = {'Inner': 0, 'Middle': 1, 'Outer': 2 } 

#title_2_sp.T[2]
#title_2_sp.T[3]

#th_min = title_2_sp.T[0].astype('float')
#th_max = title_2_sp.T[2].astype('float')
#cz_sel = title_2_sp.T[3]

th_min = n.array([ redshift_dict_min[el] for el in title_2_sp.T[2] ])
th_max = n.array([ redshift_dict_max[el] for el in title_2_sp.T[2] ])
cz_sel = n.array([ bin_dict_max[el] for el in title_2_sp.T[3] ])

#############################################
#############################################
# STEP 1, computes chi2, ndof
#############################################
#############################################

# loops over id_spec to create the 

chi2 = n.ones((len(file_list),2))*dV
ndof = n.ones((len(file_list),2))*dV
N_spectra = n.ones(len(file_list))*dV
def get_spFly(id_spec, baseName):
	"""
	Computes the chi2 and number of degrees of freedom for each model.
	"""
	path_2_spec = os.path.join( path_2_spec_dir, baseName )
	path_2_spFly = os.path.join( path_2_fly_dir, 'spFly-'+baseName )
	d = fits.open(path_2_spec)
	model = fits.open(path_2_spFly)
	converged = n.array([hdu.header['converged']=='True' for hdu in model[1:]])
	cvm = n.arange(len(model[1:]))[(converged==True)]
	if len(cvm)>=1:
		# interpolates the observations and the errors
		sel = (d[1].data['NspectraPerPixel']>0.5*n.max(d[1].data['NspectraPerPixel']))&(d[1].data['medianStack']>0)
		x_data = d[1].data['wavelength'][sel]
		y_data = d[1].data['medianStack'][sel]
		y_err = d[1].data['jackknifStackErrors'][sel]
		N_spec = d[1].data['NspectraPerPixel'][sel]

		# loops over the models to get the chi2
		for j_hdu, m in enumerate(model[1:]):
			if m.header['converged']=='True':
				x_model = m.data['wavelength']
				y_model = m.data['firefly_model']
				XMIN = n.max(n.array([n.min(x_model), n.min(x_data)]))+2.
				XMAX = n.min(n.array([n.max(x_model), n.max(x_data)]))-2.
				s_model = (x_model > XMIN-1) & (x_model < XMAX+1)
				s_data  = (x_data  > XMIN) & (x_data  < XMAX)
				itp_model = interp1d( x_model[s_model], y_model[s_model] )
				itp_data  = interp1d( x_data[s_data], y_data[s_data] )
				chi2s = abs(y_data[s_data]-itp_model(x_data[s_data]))/y_err[s_data]
				chi2[id_spec][j_hdu] =  n.sum(chi2s**2.) 
				ndof[id_spec][j_hdu] =  len(chi2s)-2.


		return 1.
	else:
		return 0.

for id_spec, baseName in enumerate(baseNames):
	get_spFly(id_spec, baseName)


#############################################
#############################################
# STEP 2
#############################################
#############################################

prihdr = fits.Header()
prihdr['file']   = os.path.basename(path_2_out_file)
prihdr['models'] = 'Maraston_2011'
#prihdr['library'] = 'MILES'
prihdr['fitter'] = 'FIREFLY'
#prihdr['author'] = 'johan comparat'
prihdr['DR'] = 16


def get_table_entry_full(hduSPM):
	"""
	Convert the results located in the headers of the spFly files into a table.
	Changes the resulting numbers to units
	"""
	if hduSPM.header['converged']=='True':
		prefix = hduSPM.header['IMF'] + "_" + hduSPM.header['MODEL'] + "_"
		table_entry = [
		  1e9*10**hduSPM.header['age_lightW'] 
		, 1e9*10**hduSPM.header['age_lightW_up_1sig']            
		, 1e9*10**hduSPM.header['age_lightW_low_1sig']           
		, 1e9*10**hduSPM.header['age_lightW_up_2sig']            
		, 1e9*10**hduSPM.header['age_lightW_low_2sig']           
		, 1e9*10**hduSPM.header['age_lightW_up_3sig']            
		, 1e9*10**hduSPM.header['age_lightW_low_3sig']           
		,     10**hduSPM.header['metallicity_lightW']            
		,     10**hduSPM.header['metallicity_lightW_up_1sig'] 
		,     10**hduSPM.header['metallicity_lightW_low_1sig']
		,     10**hduSPM.header['metallicity_lightW_up_2sig'] 
		,     10**hduSPM.header['metallicity_lightW_low_2sig']
		,     10**hduSPM.header['metallicity_lightW_up_3sig'] 
		,     10**hduSPM.header['metallicity_lightW_low_3sig']
		, 1e9*10**hduSPM.header['age_massW']                  
		, 1e9*10**hduSPM.header['age_massW_up_1sig']          
		, 1e9*10**hduSPM.header['age_massW_low_1sig']         
		, 1e9*10**hduSPM.header['age_massW_up_2sig']          
		, 1e9*10**hduSPM.header['age_massW_low_2sig']         
		, 1e9*10**hduSPM.header['age_massW_up_3sig']          
		, 1e9*10**hduSPM.header['age_massW_low_3sig']         
		,     10**hduSPM.header['metallicity_massW']          
		,     10**hduSPM.header['metallicity_massW_up_1sig']  
		,     10**hduSPM.header['metallicity_massW_low_1sig'] 
		,     10**hduSPM.header['metallicity_massW_up_2sig']    
		,     10**hduSPM.header['metallicity_massW_low_2sig']   
		,     10**hduSPM.header['metallicity_massW_up_3sig']    
		,     10**hduSPM.header['metallicity_massW_low_3sig']   
		,     10**hduSPM.header['total_mass']                   
		,     10**hduSPM.header['stellar_mass']                 
		,     10**hduSPM.header['living_stars_mass']            
		,     10**hduSPM.header['remnant_mass']                 
		,     10**hduSPM.header['remnant_mass_in_whitedwarfs']  
		,     10**hduSPM.header['remnant_mass_in_neutronstars'] 
		,     10**hduSPM.header['remnant_mass_blackholes']      
		,     10**hduSPM.header['mass_of_ejecta']     
		,     10**hduSPM.header['total_mass_up_1sig'] 
		,     10**hduSPM.header['total_mass_low_1sig']
		,     10**hduSPM.header['total_mass_up_2sig'] 
		,     10**hduSPM.header['total_mass_low_2sig']
		,     10**hduSPM.header['total_mass_up_3sig'] 
		,     10**hduSPM.header['total_mass_low_3sig']
		,         hduSPM.header['EBV']                
		,         hduSPM.header['ssp_number']         
		]
		if hduSPM.header['ssp_number'] >0 :
			ssp_num = hduSPM.header['ssp_number']
		else :
			ssp_num = 0
		#print(ssp_num)

		##print hduSPM.header
		for iii in n.arange(ssp_num):
			#header_listB = n.array([
				#' '+prefix+'total_mass_ssp_'+str(iii)                                     
				#,' '+prefix+'stellar_mass_ssp_'+str(iii)                                   
				#,' '+prefix+'living_stars_mass_ssp_'+str(iii)                              
				#,' '+prefix+'remnant_mass_ssp_'+str(iii)                                   
				#,' '+prefix+'remnant_mass_in_whitedwarfs_ssp_'+str(iii)                    
				#,' '+prefix+'remnant_mass_in_neutronstars_ssp_'+str(iii)                   
				#,' '+prefix+'remnant_mass_in_blackholes_ssp_'+str(iii)                     
				#,' '+prefix+'mass_of_ejecta_ssp_'+str(iii)                                 
				#,' '+prefix+'log_age_ssp_'+str(iii)                                        
				#,' '+prefix+'metal_ssp_'+str(iii)                                          
				#,' '+prefix+'SFR_ssp_'+str(iii)                                            
				#,' '+prefix+'weightMass_ssp_'+str(iii)                                     
				#,' '+prefix+'weightLight_ssp_'+str(iii)   
			#])
			#headerB = "".join(header_listB)
			#print(headerB)
			# values
			table_entry.append( 10**hduSPM.header['total_mass_ssp_'+str(iii)] )
			table_entry.append( 10**hduSPM.header['stellar_mass_ssp_'+str(iii)] )
			table_entry.append( 10**hduSPM.header['living_stars_mass_ssp_'+str(iii)] )
			table_entry.append( 10**hduSPM.header['remnant_mass_ssp_'+str(iii)] )
			table_entry.append( 10**hduSPM.header['remnant_mass_in_whitedwarfs_ssp_'+str(iii)] )
			table_entry.append( 10**hduSPM.header['remnant_mass_in_neutronstars_ssp_'+str(iii)] )
			table_entry.append( 10**hduSPM.header['remnant_mass_in_blackholes_ssp_'+str(iii)] )
			table_entry.append( 10**hduSPM.header['mass_of_ejecta_ssp_'+str(iii)] )
			table_entry.append( 1e9*10**hduSPM.header['log_age_ssp_'+str(iii)] )
			table_entry.append( 10**hduSPM.header['metal_ssp_'+str(iii)] )
			table_entry.append( hduSPM.header['SFR_ssp_'+str(iii)] )
			table_entry.append( hduSPM.header['weightMass_ssp_'+str(iii)] )
			table_entry.append( hduSPM.header['weightLight_ssp_'+str(iii)] )
			# concatenates headers
			#headerA += baseNamesheaderB
		
		if ssp_num<8 :
			for iii in n.arange(ssp_num, 8, 1):
				#header_listB = n.array([
					#' '+prefix+'total_mass_ssp_'+str(iii)                                     
					#,' '+prefix+'stellar_mass_ssp_'+str(iii)                                   
					#,' '+prefix+'living_stars_mass_ssp_'+str(iii)                              
					#,' '+prefix+'remnant_mass_ssp_'+str(iii)                                   
					#,' '+prefix+'remnant_mass_in_whitedwarfs_ssp_'+str(iii)                    
					#,' '+prefix+'remnant_mass_in_neutronstars_ssp_'+str(iii)                   
					#,' '+prefix+'remnant_mass_in_blackholes_ssp_'+str(iii)                     
					#,' '+prefix+'mass_of_ejecta_ssp_'+str(iii)                                 
					#,' '+prefix+'log_age_ssp_'+str(iii)                                        
					#,' '+prefix+'metal_ssp_'+str(iii)                                          
					#,' '+prefix+'SFR_ssp_'+str(iii)                                            
					#,' '+prefix+'weightMass_ssp_'+str(iii)                                     
					#,' '+prefix+'weightLight_ssp_'+str(iii)   
				#])
				#headerB = "".join(header_listB)

				table_entry.append([dV, dV, dV, dV, dV, dV, dV, dV, dV, dV, dV, dV, dV])
				#headerA += headerB

		table_entry = n.array( n.hstack((table_entry)) )
		##print table_entry.shape
		return n.hstack((table_entry))

	else:
		return n.ones(148)*dV

# step 2 : match to the created data set	

table_all = n.ones(( len(file_list), 296)) * dV
headers = ""
for index, baseName in enumerate(baseNames):
	path_2_spFly = os.path.join( path_2_fly_dir, 'spFly-'+baseName )
	d = fits.open(os.path.join(path_2_spec_dir,baseName))
	sel = (d[1].data['NspectraPerPixel'] > 0.5*n.max(d[1].data['NspectraPerPixel'])) & (d[1].data['medianStack']>0)
	#x_data = d[1].data['wavelength'][sel]
	#y_data = d[1].data['medianStack'][sel]
	#y_err = d[1].data['jackknifStackErrors'][sel]
	N_spec = d[1].data['NspectraPerPixel'][sel]
	N_spectra[index] = int(n.median(N_spec))
	hduSPM=fits.open(path_2_spFly)
	table_entry_1 = get_table_entry_full( hduSPM[1] )
	table_entry_2 = get_table_entry_full( hduSPM[2] )
	table_all[index] = n.hstack((table_entry_1, table_entry_2))

newDat = n.transpose(table_all)

headers = " Chabrier_MILES_age_lightW Chabrier_MILES_age_lightW_up_1sig Chabrier_MILES_age_lightW_low_1sig Chabrier_MILES_age_lightW_up_2sig Chabrier_MILES_age_lightW_low_2sig Chabrier_MILES_age_lightW_up_3sig Chabrier_MILES_age_lightW_low_3sig Chabrier_MILES_metallicity_lightW Chabrier_MILES_metallicity_lightW_up_1sig Chabrier_MILES_metallicity_lightW_low_1sig Chabrier_MILES_metallicity_lightW_up_2sig Chabrier_MILES_metallicity_lightW_low_2sig Chabrier_MILES_metallicity_lightW_up_3sig Chabrier_MILES_metallicity_lightW_low_3sig Chabrier_MILES_age_massW Chabrier_MILES_age_massW_up_1sig Chabrier_MILES_age_massW_low_1sig Chabrier_MILES_age_massW_up_2sig Chabrier_MILES_age_massW_low_2sig Chabrier_MILES_age_massW_up_3sig Chabrier_MILES_age_massW_low_3sig Chabrier_MILES_metallicity_massW Chabrier_MILES_metallicity_massW_up_1sig Chabrier_MILES_metallicity_massW_low_1sig Chabrier_MILES_metallicity_massW_up_2sig Chabrier_MILES_metallicity_massW_low_2sig Chabrier_MILES_metallicity_massW_up_3sig Chabrier_MILES_metallicity_massW_low_3sig Chabrier_MILES_total_mass Chabrier_MILES_stellar_mass Chabrier_MILES_living_stars_mass Chabrier_MILES_remnant_mass Chabrier_MILES_remnant_mass_in_whitedwarfs Chabrier_MILES_remnant_mass_in_neutronstars Chabrier_MILES_remnant_mass_blackholes Chabrier_MILES_mass_of_ejecta Chabrier_MILES_total_mass_up_1sig Chabrier_MILES_total_mass_low_1sig Chabrier_MILES_total_mass_up_2sig Chabrier_MILES_total_mass_low_2sig Chabrier_MILES_total_mass_up_3sig Chabrier_MILES_total_mass_low_3sig Chabrier_MILES_spm_EBV Chabrier_MILES_nComponentsSSP Chabrier_MILES_total_mass_ssp_0 Chabrier_MILES_stellar_mass_ssp_0 Chabrier_MILES_living_stars_mass_ssp_0 Chabrier_MILES_remnant_mass_ssp_0 Chabrier_MILES_remnant_mass_in_whitedwarfs_ssp_0 Chabrier_MILES_remnant_mass_in_neutronstars_ssp_0 Chabrier_MILES_remnant_mass_in_blackholes_ssp_0 Chabrier_MILES_mass_of_ejecta_ssp_0 Chabrier_MILES_log_age_ssp_0 Chabrier_MILES_metal_ssp_0 Chabrier_MILES_SFR_ssp_0 Chabrier_MILES_weightMass_ssp_0 Chabrier_MILES_weightLight_ssp_0 Chabrier_MILES_total_mass_ssp_1 Chabrier_MILES_stellar_mass_ssp_1 Chabrier_MILES_living_stars_mass_ssp_1 Chabrier_MILES_remnant_mass_ssp_1 Chabrier_MILES_remnant_mass_in_whitedwarfs_ssp_1 Chabrier_MILES_remnant_mass_in_neutronstars_ssp_1 Chabrier_MILES_remnant_mass_in_blackholes_ssp_1 Chabrier_MILES_mass_of_ejecta_ssp_1 Chabrier_MILES_log_age_ssp_1 Chabrier_MILES_metal_ssp_1 Chabrier_MILES_SFR_ssp_1 Chabrier_MILES_weightMass_ssp_1 Chabrier_MILES_weightLight_ssp_1 Chabrier_MILES_total_mass_ssp_2 Chabrier_MILES_stellar_mass_ssp_2 Chabrier_MILES_living_stars_mass_ssp_2 Chabrier_MILES_remnant_mass_ssp_2 Chabrier_MILES_remnant_mass_in_whitedwarfs_ssp_2 Chabrier_MILES_remnant_mass_in_neutronstars_ssp_2 Chabrier_MILES_remnant_mass_in_blackholes_ssp_2 Chabrier_MILES_mass_of_ejecta_ssp_2 Chabrier_MILES_log_age_ssp_2 Chabrier_MILES_metal_ssp_2 Chabrier_MILES_SFR_ssp_2 Chabrier_MILES_weightMass_ssp_2 Chabrier_MILES_weightLight_ssp_2 Chabrier_MILES_total_mass_ssp_3 Chabrier_MILES_stellar_mass_ssp_3 Chabrier_MILES_living_stars_mass_ssp_3 Chabrier_MILES_remnant_mass_ssp_3 Chabrier_MILES_remnant_mass_in_whitedwarfs_ssp_3 Chabrier_MILES_remnant_mass_in_neutronstars_ssp_3 Chabrier_MILES_remnant_mass_in_blackholes_ssp_3 Chabrier_MILES_mass_of_ejecta_ssp_3 Chabrier_MILES_log_age_ssp_3 Chabrier_MILES_metal_ssp_3 Chabrier_MILES_SFR_ssp_3 Chabrier_MILES_weightMass_ssp_3 Chabrier_MILES_weightLight_ssp_3 Chabrier_MILES_total_mass_ssp_4 Chabrier_MILES_stellar_mass_ssp_4 Chabrier_MILES_living_stars_mass_ssp_4 Chabrier_MILES_remnant_mass_ssp_4 Chabrier_MILES_remnant_mass_in_whitedwarfs_ssp_4 Chabrier_MILES_remnant_mass_in_neutronstars_ssp_4 Chabrier_MILES_remnant_mass_in_blackholes_ssp_4 Chabrier_MILES_mass_of_ejecta_ssp_4 Chabrier_MILES_log_age_ssp_4 Chabrier_MILES_metal_ssp_4 Chabrier_MILES_SFR_ssp_4 Chabrier_MILES_weightMass_ssp_4 Chabrier_MILES_weightLight_ssp_4 Chabrier_MILES_total_mass_ssp_5 Chabrier_MILES_stellar_mass_ssp_5 Chabrier_MILES_living_stars_mass_ssp_5 Chabrier_MILES_remnant_mass_ssp_5 Chabrier_MILES_remnant_mass_in_whitedwarfs_ssp_5 Chabrier_MILES_remnant_mass_in_neutronstars_ssp_5 Chabrier_MILES_remnant_mass_in_blackholes_ssp_5 Chabrier_MILES_mass_of_ejecta_ssp_5 Chabrier_MILES_log_age_ssp_5 Chabrier_MILES_metal_ssp_5 Chabrier_MILES_SFR_ssp_5 Chabrier_MILES_weightMass_ssp_5 Chabrier_MILES_weightLight_ssp_5 Chabrier_MILES_total_mass_ssp_6 Chabrier_MILES_stellar_mass_ssp_6 Chabrier_MILES_living_stars_mass_ssp_6 Chabrier_MILES_remnant_mass_ssp_6 Chabrier_MILES_remnant_mass_in_whitedwarfs_ssp_6 Chabrier_MILES_remnant_mass_in_neutronstars_ssp_6 Chabrier_MILES_remnant_mass_in_blackholes_ssp_6 Chabrier_MILES_mass_of_ejecta_ssp_6 Chabrier_MILES_log_age_ssp_6 Chabrier_MILES_metal_ssp_6 Chabrier_MILES_SFR_ssp_6 Chabrier_MILES_weightMass_ssp_6 Chabrier_MILES_weightLight_ssp_6 Chabrier_MILES_total_mass_ssp_7 Chabrier_MILES_stellar_mass_ssp_7 Chabrier_MILES_living_stars_mass_ssp_7 Chabrier_MILES_remnant_mass_ssp_7 Chabrier_MILES_remnant_mass_in_whitedwarfs_ssp_7 Chabrier_MILES_remnant_mass_in_neutronstars_ssp_7 Chabrier_MILES_remnant_mass_in_blackholes_ssp_7 Chabrier_MILES_mass_of_ejecta_ssp_7 Chabrier_MILES_log_age_ssp_7 Chabrier_MILES_metal_ssp_7 Chabrier_MILES_SFR_ssp_7 Chabrier_MILES_weightMass_ssp_7 Chabrier_MILES_weightLight_ssp_7 Chabrier_ELODIE_age_lightW Chabrier_ELODIE_age_lightW_up_1sig Chabrier_ELODIE_age_lightW_low_1sig Chabrier_ELODIE_age_lightW_up_2sig Chabrier_ELODIE_age_lightW_low_2sig Chabrier_ELODIE_age_lightW_up_3sig Chabrier_ELODIE_age_lightW_low_3sig Chabrier_ELODIE_metallicity_lightW Chabrier_ELODIE_metallicity_lightW_up_1sig Chabrier_ELODIE_metallicity_lightW_low_1sig Chabrier_ELODIE_metallicity_lightW_up_2sig Chabrier_ELODIE_metallicity_lightW_low_2sig Chabrier_ELODIE_metallicity_lightW_up_3sig Chabrier_ELODIE_metallicity_lightW_low_3sig Chabrier_ELODIE_age_massW Chabrier_ELODIE_age_massW_up_1sig Chabrier_ELODIE_age_massW_low_1sig Chabrier_ELODIE_age_massW_up_2sig Chabrier_ELODIE_age_massW_low_2sig Chabrier_ELODIE_age_massW_up_3sig Chabrier_ELODIE_age_massW_low_3sig Chabrier_ELODIE_metallicity_massW Chabrier_ELODIE_metallicity_massW_up_1sig Chabrier_ELODIE_metallicity_massW_low_1sig Chabrier_ELODIE_metallicity_massW_up_2sig Chabrier_ELODIE_metallicity_massW_low_2sig Chabrier_ELODIE_metallicity_massW_up_3sig Chabrier_ELODIE_metallicity_massW_low_3sig Chabrier_ELODIE_total_mass Chabrier_ELODIE_stellar_mass Chabrier_ELODIE_living_stars_mass Chabrier_ELODIE_remnant_mass Chabrier_ELODIE_remnant_mass_in_whitedwarfs Chabrier_ELODIE_remnant_mass_in_neutronstars Chabrier_ELODIE_remnant_mass_blackholes Chabrier_ELODIE_mass_of_ejecta Chabrier_ELODIE_total_mass_up_1sig Chabrier_ELODIE_total_mass_low_1sig Chabrier_ELODIE_total_mass_up_2sig Chabrier_ELODIE_total_mass_low_2sig Chabrier_ELODIE_total_mass_up_3sig Chabrier_ELODIE_total_mass_low_3sig Chabrier_ELODIE_spm_EBV Chabrier_ELODIE_nComponentsSSP Chabrier_ELODIE_total_mass_ssp_0 Chabrier_ELODIE_stellar_mass_ssp_0 Chabrier_ELODIE_living_stars_mass_ssp_0 Chabrier_ELODIE_remnant_mass_ssp_0 Chabrier_ELODIE_remnant_mass_in_whitedwarfs_ssp_0 Chabrier_ELODIE_remnant_mass_in_neutronstars_ssp_0 Chabrier_ELODIE_remnant_mass_in_blackholes_ssp_0 Chabrier_ELODIE_mass_of_ejecta_ssp_0 Chabrier_ELODIE_log_age_ssp_0 Chabrier_ELODIE_metal_ssp_0 Chabrier_ELODIE_SFR_ssp_0 Chabrier_ELODIE_weightMass_ssp_0 Chabrier_ELODIE_weightLight_ssp_0 Chabrier_ELODIE_total_mass_ssp_1 Chabrier_ELODIE_stellar_mass_ssp_1 Chabrier_ELODIE_living_stars_mass_ssp_1 Chabrier_ELODIE_remnant_mass_ssp_1 Chabrier_ELODIE_remnant_mass_in_whitedwarfs_ssp_1 Chabrier_ELODIE_remnant_mass_in_neutronstars_ssp_1 Chabrier_ELODIE_remnant_mass_in_blackholes_ssp_1 Chabrier_ELODIE_mass_of_ejecta_ssp_1 Chabrier_ELODIE_log_age_ssp_1 Chabrier_ELODIE_metal_ssp_1 Chabrier_ELODIE_SFR_ssp_1 Chabrier_ELODIE_weightMass_ssp_1 Chabrier_ELODIE_weightLight_ssp_1 Chabrier_ELODIE_total_mass_ssp_2 Chabrier_ELODIE_stellar_mass_ssp_2 Chabrier_ELODIE_living_stars_mass_ssp_2 Chabrier_ELODIE_remnant_mass_ssp_2 Chabrier_ELODIE_remnant_mass_in_whitedwarfs_ssp_2 Chabrier_ELODIE_remnant_mass_in_neutronstars_ssp_2 Chabrier_ELODIE_remnant_mass_in_blackholes_ssp_2 Chabrier_ELODIE_mass_of_ejecta_ssp_2 Chabrier_ELODIE_log_age_ssp_2 Chabrier_ELODIE_metal_ssp_2 Chabrier_ELODIE_SFR_ssp_2 Chabrier_ELODIE_weightMass_ssp_2 Chabrier_ELODIE_weightLight_ssp_2 Chabrier_ELODIE_total_mass_ssp_3 Chabrier_ELODIE_stellar_mass_ssp_3 Chabrier_ELODIE_living_stars_mass_ssp_3 Chabrier_ELODIE_remnant_mass_ssp_3 Chabrier_ELODIE_remnant_mass_in_whitedwarfs_ssp_3 Chabrier_ELODIE_remnant_mass_in_neutronstars_ssp_3 Chabrier_ELODIE_remnant_mass_in_blackholes_ssp_3 Chabrier_ELODIE_mass_of_ejecta_ssp_3 Chabrier_ELODIE_log_age_ssp_3 Chabrier_ELODIE_metal_ssp_3 Chabrier_ELODIE_SFR_ssp_3 Chabrier_ELODIE_weightMass_ssp_3 Chabrier_ELODIE_weightLight_ssp_3 Chabrier_ELODIE_total_mass_ssp_4 Chabrier_ELODIE_stellar_mass_ssp_4 Chabrier_ELODIE_living_stars_mass_ssp_4 Chabrier_ELODIE_remnant_mass_ssp_4 Chabrier_ELODIE_remnant_mass_in_whitedwarfs_ssp_4 Chabrier_ELODIE_remnant_mass_in_neutronstars_ssp_4 Chabrier_ELODIE_remnant_mass_in_blackholes_ssp_4 Chabrier_ELODIE_mass_of_ejecta_ssp_4 Chabrier_ELODIE_log_age_ssp_4 Chabrier_ELODIE_metal_ssp_4 Chabrier_ELODIE_SFR_ssp_4 Chabrier_ELODIE_weightMass_ssp_4 Chabrier_ELODIE_weightLight_ssp_4 Chabrier_ELODIE_total_mass_ssp_5 Chabrier_ELODIE_stellar_mass_ssp_5 Chabrier_ELODIE_living_stars_mass_ssp_5 Chabrier_ELODIE_remnant_mass_ssp_5 Chabrier_ELODIE_remnant_mass_in_whitedwarfs_ssp_5 Chabrier_ELODIE_remnant_mass_in_neutronstars_ssp_5 Chabrier_ELODIE_remnant_mass_in_blackholes_ssp_5 Chabrier_ELODIE_mass_of_ejecta_ssp_5 Chabrier_ELODIE_log_age_ssp_5 Chabrier_ELODIE_metal_ssp_5 Chabrier_ELODIE_SFR_ssp_5 Chabrier_ELODIE_weightMass_ssp_5 Chabrier_ELODIE_weightLight_ssp_5 Chabrier_ELODIE_total_mass_ssp_6 Chabrier_ELODIE_stellar_mass_ssp_6 Chabrier_ELODIE_living_stars_mass_ssp_6 Chabrier_ELODIE_remnant_mass_ssp_6 Chabrier_ELODIE_remnant_mass_in_whitedwarfs_ssp_6 Chabrier_ELODIE_remnant_mass_in_neutronstars_ssp_6 Chabrier_ELODIE_remnant_mass_in_blackholes_ssp_6 Chabrier_ELODIE_mass_of_ejecta_ssp_6 Chabrier_ELODIE_log_age_ssp_6 Chabrier_ELODIE_metal_ssp_6 Chabrier_ELODIE_SFR_ssp_6 Chabrier_ELODIE_weightMass_ssp_6 Chabrier_ELODIE_weightLight_ssp_6 Chabrier_ELODIE_total_mass_ssp_7 Chabrier_ELODIE_stellar_mass_ssp_7 Chabrier_ELODIE_living_stars_mass_ssp_7 Chabrier_ELODIE_remnant_mass_ssp_7 Chabrier_ELODIE_remnant_mass_in_whitedwarfs_ssp_7 Chabrier_ELODIE_remnant_mass_in_neutronstars_ssp_7 Chabrier_ELODIE_remnant_mass_in_blackholes_ssp_7 Chabrier_ELODIE_mass_of_ejecta_ssp_7 Chabrier_ELODIE_log_age_ssp_7 Chabrier_ELODIE_metal_ssp_7 Chabrier_ELODIE_SFR_ssp_7 Chabrier_ELODIE_weightMass_ssp_7 Chabrier_ELODIE_weightLight_ssp_7"



































































t = Table()

for data_array, head in zip(newDat, headers.split()):
	t.add_column(Column(name=head, format='D', data=data_array))

for id_col, (col_chi2, col_ndof) in enumerate(zip(n.transpose(chi2), n.transpose(ndof))):
	t.add_column(Column(name=hdu_header_prefix[id_col]+"chi2", format='D', data=col_chi2))
	t.add_column(Column(name=hdu_header_prefix[id_col]+"ndof", format='D', data=col_ndof))

t.add_column(Column(name="z_min", format='D', data=th_min))
t.add_column(Column(name="z_max", format='D', data=th_max))
t.add_column(Column(name="radial_bin", format='D', data=cz_sel))
t.add_column(Column(name="Inner", format='L', data=(cz_sel==0)))
t.add_column(Column(name="Middle", format='L', data=(cz_sel==1)))
t.add_column(Column(name="Outer", format='L', data=(cz_sel==2)))

#t.add_column(Column(name="min_log10_theta_over_theta200c", format='D', data=th_min))
#t.add_column(Column(name="max_log10_theta_over_theta200c", format='D', data=th_max))
#t.add_column(Column(name="log10_theta_over_theta200c", format='D', data=0.5*(th_min+th_max)))
#t.add_column(Column(name="CZ_selection", format='12A', data=cz_sel))
#t.add_column(Column(name="allCZ", format='L', data=cz_sel=='allCZ'))
#t.add_column(Column(name="highCZ", format='L', data=cz_sel=='highCZ'))
#t.add_column(Column(name="lowCZ", format='L', data=cz_sel=='lowCZ'))
t.add_column(Column(name="N_spectra", format='K', data=N_spectra))

t.write(path_2_out_file, overwrite=True)

sys.exit()
new_cols = fits.ColDefs(all_cols)
tt = Table(new_cols)
sys.exit()
tbhdu = fits.BinTableHDU.from_columns( new_cols)

#prihdr['ageMin'] = 0  
#prihdr['ageMax'] = 20 
#prihdr['Zmin']   = 0.001 
#prihdr['Zmax']   = 10  

prihdu = fits.PrimaryHDU(header=prihdr)

hdu = fits.HDUList([prihdu, tbhdu])

if os.path.isfile(path_2_out_file):
    os.remove(path_2_out_file)

hdu.writeto(path_2_out_file)


