import os

topdir = '/data36s/comparat/SDSS/ELG/clustering/'

def write_param(paramfile, data, random, result):
	f=open(paramfile, 'w')
	f.write('data_filename= '+data+' \n')     
	f.write('random_filename= '+random+' \n')
	f.write('input_format= 2  \n')
	f.write('output_filename= '+result+' \n')
	f.write('corr_type= monopole \n')
	f.write('omega_M= 0.307 \n')
	f.write('omega_L= 0.693 \n')
	f.write('w= -1          \n')
	f.write('log_bin= 0     \n')
	f.write('dim1_max= 200. \n')
	#f.write('dim1_min_logbin= 0.2 \n')
	f.write('dim1_nbin= 50  \n')
	f.write('dim2_max= 200. \n')
	f.write('dim2_nbin= 50  \n')
	f.write('dim3_min= 0.7 \n')
	f.write('dim3_max= 1.1  \n')
	f.write('dim3_nbin= 1   \n')
	f.close()


write_param("ebosselg.JC28Nov2017_data.param"                   , topdir+"ebosselg.JC28Nov2017_data"                        , topdir+"ebosselg.JC28Nov2017_rand"                        , topdir+"ebosselg.JC28Nov2017_linB4_2pcf"                         )
write_param("ebosselg.JC28Nov2017_data.mask_O2O3Hb_th_1.5.param", topdir+"ebosselg.JC28Nov2017_data.mask_O2O3Hb_th_1.5.txt" , topdir+"ebosselg.JC28Nov2017_rand.mask_O2O3Hb_th_1.5.txt" , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.mask_O2O3Hb_th_1.5.txt"  )
write_param("ebosselg.JC28Nov2017_data.mask_O2O3Hb_th_1.param"  , topdir+"ebosselg.JC28Nov2017_data.mask_O2O3Hb_th_1.txt"   , topdir+"ebosselg.JC28Nov2017_rand.mask_O2O3Hb_th_1.txt"   , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.mask_O2O3Hb_th_1.txt"    )
write_param("ebosselg.JC28Nov2017_data.mask_O2O3Hb_th_2.param"  , topdir+"ebosselg.JC28Nov2017_data.mask_O2O3Hb_th_2.txt"   , topdir+"ebosselg.JC28Nov2017_rand.mask_O2O3Hb_th_2.txt"   , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.mask_O2O3Hb_th_2.txt"    )
write_param("ebosselg.JC28Nov2017_data.mask_O2O3_th_1.5.param"  , topdir+"ebosselg.JC28Nov2017_data.mask_O2O3_th_1.5.txt"   , topdir+"ebosselg.JC28Nov2017_rand.mask_O2O3_th_1.5.txt"   , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.mask_O2O3_th_1.5.txt"    )
write_param("ebosselg.JC28Nov2017_data.mask_O2O3_th_1.param"    , topdir+"ebosselg.JC28Nov2017_data.mask_O2O3_th_1.txt"     , topdir+"ebosselg.JC28Nov2017_rand.mask_O2O3_th_1.txt"     , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.mask_O2O3_th_1.txt"      )
write_param("ebosselg.JC28Nov2017_data.mask_O2O3_th_2.param"    , topdir+"ebosselg.JC28Nov2017_data.mask_O2O3_th_2.txt"     , topdir+"ebosselg.JC28Nov2017_rand.mask_O2O3_th_2.txt"     , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.mask_O2O3_th_2.txt"      )
write_param("ebosselg.JC28Nov2017_data.mask_O2_th_1.5.param"    , topdir+"ebosselg.JC28Nov2017_data.mask_O2_th_1.5.txt"     , topdir+"ebosselg.JC28Nov2017_rand.mask_O2_th_1.5.txt"     , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.mask_O2_th_1.5.txt"      )
write_param("ebosselg.JC28Nov2017_data.mask_O2_th_1.param"      , topdir+"ebosselg.JC28Nov2017_data.mask_O2_th_1.txt"       , topdir+"ebosselg.JC28Nov2017_rand.mask_O2_th_1.txt"       , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.mask_O2_th_1.txt"        )
write_param("ebosselg.JC28Nov2017_data.mask_O2_th_2.param"      , topdir+"ebosselg.JC28Nov2017_data.mask_O2_th_2.txt"       , topdir+"ebosselg.JC28Nov2017_rand.mask_O2_th_2.txt"       , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.mask_O2_th_2.txt"        )
write_param("ebosselg.JC28Nov2017_data.O2O3Hb_th_3.param"       , topdir+"ebosselg.JC28Nov2017_data.O2O3Hb_th_3.txt"        , topdir+"ebosselg.JC28Nov2017_rand.O2O3Hb_th_3.txt"        , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.O2O3Hb_th_3.txt"         )
write_param("ebosselg.JC28Nov2017_data.O2O3_th_3.param"         , topdir+"ebosselg.JC28Nov2017_data.O2O3_th_3.txt"          , topdir+"ebosselg.JC28Nov2017_rand.O2O3_th_3.txt"          , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.O2O3_th_3.txt"           )
write_param("ebosselg.JC28Nov2017_data.O2_th_3.param"           , topdir+"ebosselg.JC28Nov2017_data.O2_th_3.txt"            , topdir+"ebosselg.JC28Nov2017_rand.O2_th_3.txt"            , topdir+"ebosselg.JC28Nov2017_linB4_2pcf.O2_th_3.txt"             )















