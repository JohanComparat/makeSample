"""
Very important class !

Setup the snapshot replication in x, y, z.

Calculates automatically all the rotations and translations.

All other scripts depend on this !

"""
class lc_setup_class(object):
  pass

from scipy.interpolate import interp1d
import numpy as n
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)

z_array = n.arange(0, 7., 0.005)
dc_to_z = interp1d(cosmoMD.comoving_distance(z_array), z_array)

from scipy.integrate import quad

def compute_area(th_min, th_max, phi_min, phi_max):
	area =quad(n.sin, th_min*n.pi/180., th_max*n.pi/180.)[0] * (phi_max*n.pi/180. - phi_min*n.pi/180.) * (180./n.pi)**2.
	print(area, 'deg2')
	return area

#z_min = dc_to_z(d_min)#.value
#z_max = dc_to_z(d_max)#.value
"""
print("dimensions of the light cone")

C2 = lc_setup()
C4 = lc_setup()
C6 = lc_setup()
L_box = 1000./0.6777

C2.name = 'C2'
C2.remap_los = 'x'
C2.x_obs, C2.y_obs, C2.z_obs = 0., 1.0000/2.*L_box, 0.7071/2.*L_box
C2.strech_factor_los = 1.4142 
C2.strech_factor_y = 1.0
C2.strech_factor_z = 0.7071

C2.d_min = 0. 
C2.d_max = C2.d_min + C2.strech_factor_los * L_box

C2.theta_max = n.arctan2(C2.strech_factor_y*L_box,C2.d_max)*180./n.pi
C2.phi_max = n.arctan2(C2.strech_factor_z*L_box,C2.d_max)*180./n.pi

C2.z_start = dc_to_z(C2.d_min)
C2.z_reach = dc_to_z(C2.d_max)
print(C2.name, C2.z_start,"z<",C2.z_reach,"|theta [deg]|<",C2.theta_max, "|phi [deg]|<", C2.phi_max)

#C2.z_start = dc_to_z(C2.d_min)
#C2.z_reach = dc_to_z(C2.d_max)
#C2.dec_max = (C2.strech_factor_y*L_box*u.Mpc/cosmoMD.kpc_comoving_per_arcmin(C2.z_reach).to(u.Mpc/u.degree)/2.).value
#C2.ra_max =  (C2.strech_factor_z*L_box*u.Mpc/cosmoMD.kpc_comoving_per_arcmin(C2.z_reach).to(u.Mpc/u.degree)/2.).value
#C2.area = C2.ra_max*C2.dec_max*2. * 2
#print(C2.name, C2.z_start,"z<",C2.z_reach,"|ra [deg]|<",C2.ra_max, "|dec [deg]|<", C2.dec_max)

C4.name ='C4'
C4.remap_los = 'y'
C4.x_obs, C4.y_obs, C4.z_obs = -1.4142*L_box, 1.4142/2.*L_box,  0.5774/2.*L_box

C4.strech_factor_los = 1.2247
C4.strech_factor_y = 1.4142
C4.strech_factor_z = 0.5774

C4.d_min = n.min(n.array([C2.d_max*n.cos(C2.theta_max*n.pi/180.), C2.d_max*n.cos(C2.phi_max*n.pi/180.)]))
C4.d_max = C4.d_min + C4.strech_factor_los * L_box 

C4.theta_max = n.arctan2(C4.strech_factor_y*L_box,C4.d_max)*180./n.pi
C4.phi_max = n.arctan2(C4.strech_factor_z*L_box,C4.d_max)*180./n.pi

C4.z_start = dc_to_z(C4.d_min)
C4.z_reach = dc_to_z(C4.d_max)
print(C4.name, C4.z_start,"z<",C4.z_reach,"|theta [deg]|<",C4.theta_max, "|phi [deg]|<", C4.phi_max)


C6.name='C6'
#1.2112310321535573 z< 2.929816027613037 |ra [deg]|< 4.634390951003509 |dec [deg]|< 5.35140745509032
C6.remap_los = 'x'
C6.x_obs, C6.y_obs, C6.z_obs = -(1.4142+1.2247)*L_box, 0.8165/2.*L_box, 0.7071/2.*L_box

C6.strech_factor_los = 1.7321
C6.strech_factor_y = 0.8165
C6.strech_factor_z = 0.7071

C6.d_min = n.min(n.array([C4.d_max*n.cos(C4.theta_max*n.pi/180.), C4.d_max*n.cos(C4.phi_max*n.pi/180.)]))
C6.d_max = C6.d_min + C6.strech_factor_los * L_box 

C6.theta_max = n.arctan2(C6.strech_factor_y*L_box,C6.d_max)*180./n.pi
C6.phi_max = n.arctan2(C6.strech_factor_z*L_box,C6.d_max)*180./n.pi

C6.z_start = dc_to_z(C6.d_min)
C6.z_reach = dc_to_z(C6.d_max)
print(C6.name, C6.z_start,"z<",C6.z_reach,"|theta [deg]|<",C6.theta_max, "|phi [deg]|<", C6.phi_max)

n.min(n.array([C6.d_max*n.cos(C6.theta_max*n.pi/180.), C6.d_max*n.cos(C6.phi_max*n.pi/180.)]))



C6.d_min_th = (1.4142+1.2247)*L_box
C6.d_max_th = (1.4142+1.2247+1.7321)*L_box 
C6.z_start_th = dc_to_z(C6.d_min_th)
C6.z_reach_th = dc_to_z(C6.d_max_th)
C6.dec_max_th = (C6.strech_factor_y*L_box*u.Mpc/cosmoMD.kpc_comoving_per_arcmin(C6.z_reach_th).to(u.Mpc/u.degree)/2.).value
C6.ra_max_th =  (C6.strech_factor_z*L_box*u.Mpc/cosmoMD.kpc_comoving_per_arcmin(C6.z_reach_th).to(u.Mpc/u.degree)/2.).value
#print(C6.name, C6.z_start_th,"z<",C6.z_reach_th,"|ra [deg]|<",C6.ra_max_th, "|dec [deg]|<", C6.dec_max_th)

C6.Dproj_ra = n.cos(ra_max*n.pi/180.) * d_min_th
C6.Dproj_dec = n.cos(dec_max*n.pi/180.) * d_min_th
C6.d_min = n.min([Dproj_ra, Dproj_dec])
C6.d_max = d_max_th - (d_min_th-d_min)


C6.Dproj_ra  = n.cos(C6.ra_max_th*n.pi/180.) *  C6.d_min_th
C6.Dproj_dec = n.cos(C6.dec_max_th*n.pi/180.) * C6.d_min_th
C6.d_min = n.min([C6.Dproj_ra, C6.Dproj_dec])
C6.d_max = C6.d_max_th - (C6.d_min_th-C6.d_min)

C6.z_start = dc_to_z(C6.d_min)
C6.z_reach = dc_to_z(C6.d_max)
C6.dec_max = (C6.strech_factor_y*L_box*u.Mpc/cosmoMD.kpc_comoving_per_arcmin(C6.z_reach).to(u.Mpc/u.degree)/2.).value
C6.ra_max =  (C6.strech_factor_z*L_box*u.Mpc/cosmoMD.kpc_comoving_per_arcmin(C6.z_reach).to(u.Mpc/u.degree)/2.).value
C6.area = C6.ra_max*C6.dec_max*2. * 2
print(C6.name, C6.z_start,"z<",C6.z_reach,"|ra [deg]|<",C6.ra_max, "|dec [deg]|<", C6.dec_max)


"""
def get_ra_dec_bounds(setup):
	setup.theta_bounds = n.array([n.arccos( setup.z_max / setup.d_max ) * 180/n.pi, n.arccos( setup.z_min / setup.d_max ) * 180/n.pi])
	setup.phi_bounds = n.array([n.arctan2( setup.y_max , setup.d_max ) * 180/n.pi, n.arctan2( setup.y_min , setup.d_max ) * 180/n.pi])
	print(setup.theta_bounds, setup.phi_bounds)

C1 = lc_setup_class()
L_box = 1000./0.6777

C1.name = 'C1'
C1.remap_los = 'x'
C1.x_obs, C1.y_obs, C1.z_obs = 0., -0.5*L_box, -0.5*L_box
C1.strech_factor_los = 1.0
C1.strech_factor_y = 1.0
C1.strech_factor_z = 1.0

C1.x_max = L_box + C1.x_obs
C1.y_max = L_box + C1.y_obs
C1.z_max = L_box + C1.z_obs
C1.x_min = C1.x_obs
C1.y_min = C1.y_obs
C1.z_min = C1.z_obs

C1.d_min = 0. 
C1.d_start = 0. 
C1.d_max = C1.d_min + C1.strech_factor_los * L_box

C1.theta_max = n.arctan2(C1.strech_factor_y*L_box/2.,C1.d_max)*180./n.pi
C1.phi_max = n.arctan2(C1.strech_factor_z*L_box/2.,C1.d_max)*180./n.pi

C1.z_start = dc_to_z(C1.d_min)
C1.z_reach = dc_to_z(C1.d_max)
print(C1.d_min/300., C1.d_min, C1.d_max)
print(C1.name, C1.z_start,"z<",C1.z_reach,"|theta [deg]|<",C1.theta_max, "|phi [deg]|<", C1.phi_max)

get_ra_dec_bounds(C1)
C1.area = compute_area(C1.theta_bounds.min(), C1.theta_bounds.max(), C1.phi_bounds.min(), C1.phi_bounds.max())

def create_np1(C2, name):
  C4 = lc_setup_class()
  C4.name = name
  C4.remap_los = 'y'
  #
  C4.strech_factor_los = 1.0
  C4.strech_factor_y   = 1.0
  C4.strech_factor_z   = 1.0
  #
  C4.d_start = C2.d_max
  C4.d_min = n.min(n.array([C2.d_max*n.cos(C2.theta_max*n.pi/180.), C2.d_max*n.cos(C2.phi_max*n.pi/180.)]))
  C4.d_max = C4.d_min + C4.strech_factor_los * L_box 
  #
  C4.theta_max = n.arctan2(C4.strech_factor_y*L_box/2.,C4.d_max)*180./n.pi
  C4.phi_max = n.arctan2(C4.strech_factor_z*L_box/2.,C4.d_max)*180./n.pi
  #
  C4.z_start = dc_to_z(C4.d_min)
  C4.z_reach = dc_to_z(C4.d_max)
  print(C4.d_min, C4.d_max)
  print(C4.z_start,"z<",C4.z_reach,"|theta [deg]|<",C4.theta_max, "|phi [deg]|<", C4.phi_max)
  return C4



C2 = create_np1(C1, 'C2')
C2.x_obs, C2.y_obs, C2.z_obs = C2.d_min, -0.5*L_box, -0.5*L_box
C2.remap_los = 'x'
C2.x_max = L_box + C2.x_obs
C2.y_max = L_box + C2.y_obs
C2.z_max = L_box + C2.z_obs
C2.x_min = C2.x_obs
C2.y_min = C2.y_obs
C2.z_min = C2.z_obs
get_ra_dec_bounds(C2)
C2.area = compute_area(C2.theta_bounds.min(), C2.theta_bounds.max(), C2.phi_bounds.min(), C2.phi_bounds.max())

C3 = create_np1(C2, 'C3')
C3.x_obs, C3.y_obs, C3.z_obs = C3.d_min, -0.5*L_box, -0.5*L_box
C3.remap_los = 'x'
C3.x_max = L_box + C3.x_obs
C3.y_max = L_box + C3.y_obs
C3.z_max = L_box + C3.z_obs
C3.x_min = C3.x_obs
C3.y_min = C3.y_obs
C3.z_min = C3.z_obs
get_ra_dec_bounds(C3)
C3.area = compute_area(C3.theta_bounds.min(), C3.theta_bounds.max(), C3.phi_bounds.min(), C3.phi_bounds.max())

C4 = create_np1(C3, 'C4')
C4.x_obs, C4.y_obs, C4.z_obs = C4.d_min, -0.5*L_box, -0.5*L_box
C4.remap_los = 'x'
C4.x_max = L_box + C4.x_obs
C4.y_max = L_box + C4.y_obs
C4.z_max = L_box + C4.z_obs
C4.x_min = C4.x_obs
C4.y_min = C4.y_obs
C4.z_min = C4.z_obs
get_ra_dec_bounds(C4)
C4.area = compute_area(C4.theta_bounds.min(), C4.theta_bounds.max(), C4.phi_bounds.min(), C4.phi_bounds.max())

C5 = create_np1(C4, 'C5')
C5.x_obs, C5.y_obs, C5.z_obs = C5.d_min, -0.5*L_box, -0.5*L_box
C5.remap_los = 'x'
C5.x_max = L_box + C5.x_obs
C5.y_max = L_box + C5.y_obs
C5.z_max = L_box + C5.z_obs
C5.x_min = C5.x_obs
C5.y_min = C5.y_obs
C5.z_min = C5.z_obs
get_ra_dec_bounds(C5)
C5.area = compute_area(C5.theta_bounds.min(), C5.theta_bounds.max(), C5.phi_bounds.min(), C5.phi_bounds.max())

C6 = create_np1(C5, 'C6')
C6.x_obs, C6.y_obs, C6.z_obs = C6.d_min, -0.5*L_box, -0.5*L_box
C6.remap_los = 'x'
C6.x_max = L_box + C6.x_obs
C6.y_max = L_box + C6.y_obs
C6.z_max = L_box + C6.z_obs
C6.x_min = C6.x_obs
C6.y_min = C6.y_obs
C6.z_min = C6.z_obs
get_ra_dec_bounds(C6)
C6.area = compute_area(C6.theta_bounds.min(), C6.theta_bounds.max(), C6.phi_bounds.min(), C6.phi_bounds.max())


cones = {'C1': C1, 'C2': C2, 'C3': C3, 'C4': C4, 'C5': C5, 'C6': C6 }


