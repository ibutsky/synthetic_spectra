from astropy.io import fits
import numpy as np
import glob
import os
import sys
#import pynbody
import yt
from yt import YTArray
from yt.units import YTQuantity

def rotate(v, phi, theta):
    y_rot_mat = [
        [np.cos(theta),  0, np.sin(theta)],
        [0,              1, 0],
        [-np.sin(theta), 0, np.cos(theta)]]
    z_rot_mat = [
        [np.cos(phi), -np.sin(phi), 0],
        [np.sin(phi), np.cos(phi),  0],
        [0,           0,            1]]
    rot_mat = np.matmul(z_rot_mat, y_rot_mat)
    return np.matmul(rot_mat, v)

def assert_within(value, target, error):
    assert(value < target + error)
    assert(value > target - error)

def generate_random_ray_coordinates(center, rmin, rmax, ray_len):
    # center: the coordinates of the center of the simulation
    # rmin, rmax: the minimum and maximum impact parameter values of the generated ray
    # ray_len: the length of the ray, in kpc

    # random theta and phi
    theta = np.pi*np.random.random()
    phi = 2.*np.pi*np.random.random()

    # calculate r, random impact parameter in the range (rmin, rmax)
    r = (rmax-rmin)*np.random.random() + rmin

    # [0, 0, 1] is the initial dummy direction of the vector
    # named the "chosen basis" below

    origin_position = r*rotate([0,0,1], phi, theta)

    # Choose an angle in the tanget plane
    psi = 2 * np.pi * np.random.random()
    # Create unit vector along this angle in the basis of \hat{z}' = (x,y,z)
    # Call it the "chosen basis"
    z_ray_end = ray_len / 2 * np.array([np.cos(psi), np.sin(psi), 0])
    z_ray_start = ray_len / 2 * np.array([np.cos(psi + np.pi), np.sin(psi + np.pi), 0])
    assert_within(np.dot(z_ray_end, z_ray_start), -(ray_len/2)**2, 0.00001)
    # Rotate from chosen basis into standard basis        

    origin_ray_end = rotate(z_ray_end, phi, theta)
    origin_ray_start = rotate(z_ray_start, phi, theta)
    assert_within(np.dot(origin_ray_end, origin_position), 0, 0.00001)
    assert_within(np.dot(origin_ray_start, origin_position), 0, 0.00001)

    # Translate vector from the origin
    position = origin_position + center
    ray_start = origin_ray_start + position
    ray_end = origin_ray_end + position
    assert_within(np.dot(position-ray_start, position-ray_end), -(ray_len/2)**2, 0.00001)

    return r, ray_start, ray_end

def get_rockstar_data(rstar_fn, halo_id):
    """
    Use ytree to get all of the halo centroids, virial radii, and redshift info; store in a dict
    """
    # load up dataset and get appropriate TreeNode
    import ytree 

    a = ytree.load(rstar_fn)
    t = a[a["Orig_halo_ID"] == halo_id][0]

    redshift_arr = t['prog', 'redshift']
    x_arr = t['prog', 'x'].in_units('unitary')
    y_arr = t['prog', 'y'].in_units('unitary')
    z_arr = t['prog', 'z'].in_units('unitary')
    vx_arr = t['prog', 'vx'].in_units('km/s')
    vy_arr = t['prog', 'vy'].in_units('km/s')
    vz_arr = t['prog', 'vz'].in_units('km/s')
    rvir_arr = t['prog', 'Rvir'].convert_to_units('kpc') 

    return {'redshift_arr':redshift_arr, 'x_arr':x_arr, 'y_arr':y_arr, 'z_arr':z_arr, 'vx_arr':vx_arr, 'vy_arr':vy_arr, 'vz_arr':vz_arr, 'rvir_arr':rvir_arr}

def read_rockstar_info(rockstar_data, ds):
    """
    Interpolate halo center from rockstar merger tree
    """
    redshift = ds.current_redshift
    redshift_arr = rockstar_data['redshift_arr']
    x = np.interp(redshift, redshift_arr, rockstar_data['x_arr'].in_units('unitary'))
    y = np.interp(redshift, redshift_arr, rockstar_data['y_arr'].in_units('unitary'))
    z = np.interp(redshift, redshift_arr, rockstar_data['z_arr'].in_units('unitary'))
    vx = np.interp(redshift, redshift_arr, rockstar_data['vx_arr'].in_units('km/s'))
    vy = np.interp(redshift, redshift_arr, rockstar_data['vy_arr'].in_units('km/s'))
    vz = np.interp(redshift, redshift_arr, rockstar_data['vz_arr'].in_units('km/s'))
    rvir = np.interp(redshift, redshift_arr, rockstar_data['rvir_arr'])

    pos_arr = ds.arr([x,y,z], 'unitary')
    vel_arr = ds.arr([vx,vy,vz], 'km/s')
    rvir = ds.quan(rvir, 'kpccm').in_units('kpc')
    return pos_arr, vel_arr, rvir
    
if __name__ == '__main__':

    snapshot_fn = sys.argv[1]
    rockstar_fn = sys.argv[2]
    halo_id = int(sys.argv[3])
    rockstar_data = get_rockstar_data(rockstar_fn, halo_id)

    ds = yt.load(snapshot_fn)
    c, bv, rvir = read_rockstar_info(rockstar_data, ds)
    print(c)
    print(bv)
    print(rvir)


def stitch_g130m_g160m_spectra(fn_g130, fn_g160, fn_combined):
    temp130 = fits.open(fn_g130)
    temp160 = fits.open(fn_g160)
    data = []
    for item in temp130[1].columns:
        if item.name in ['wavelength', 'flux', 'flux_error']:
            # this iterates through each of the arrays in the extension by grabbing names from the header 
            data1     = temp130[1].data[item.name] # get data from file 1
            data2     = temp160[1].data[item.name]
            new_array = np.concatenate([data1, data2]) # make the new array composed of both parts 
            data.append(new_array)

    c1 = fits.Column(name = 'WAVELENGTH', array = data[0], format = 'E')
    c2 = fits.Column(name = 'FLUX', array = data[1], format = 'E')
    c3 = fits.Column(name = 'FLUXERR', array = data[2], format = 'E')
    fits_combined = fits.BinTableHDU.from_columns([c1, c2, c3])
    fits_combined.writeto(fn_combined, overwrite = True)


def load_simulation_properties(model, output, laptop = False):
    if model == 'tempest':
        ds = yt.load('/mnt/c/scratch/sciteam/chummels/Tempest10/DD%04d/DD%04d'%(output, output))
        if output == 524:
            rockstar_fn = '/projects/eot/bafa/tempest/tree_27.dat'
            halo_id = 27
            rockstar_data = get_rockstar_data(rockstar_fn, halo_id)
            gcenter, bulk_velocity, rvir = read_rockstar_info(rockstar_data, ds)
        else:
            print('WARNING: No rockstar file for output %i'%(output))
            ad = ds.all_data()
            # TODO 
    elif model == 'P0':
        # calculate center of mass and bulk velocity using pynbody
        if laptop:
            fn  = '~/Work/galaxy/P0/P0.003195'
            ds = yt.load(fn)
            gcenter = YTArray([-1.693207e+04, -1.201068e+04, 5.303337e+03], 'kpc')
            bulk_velocity = YTArray([72.78, -248.83, -63.49], 'km/s')
        else:
            import pynbody
            fn = '/nobackup/ibutsky/tmp/pioneer.%06d'%(output)
            ds = yt.load(fn)
            pynbody_file = '/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.%06d'%(output)
            s = pynbody.load(pynbody_file)
            s.physical_units()
            gcenter = YTArray(pynbody.analysis.halo.center_of_mass(s.s), 'kpc')
            gcenter = YTArray([-1.693207e+04, -1.201068e4, 5.303337e3], 'kpc')
            print(gcenter)
            bulk_velocity = YTArray(pynbody.analysis.halo.center_of_mass_velocity(s.g), 'km/s')
    
    return ds, gcenter, bulk_velocity
            

def get_next_ray_id(model, redshift, spectrum_directory = '.'):
    fn = '%s/%s_z%0.2f_ray_data.dat'%(spectrum_directory, model, redshift)
    if os.path.isfile(fn):
        ray_id_list = np.loadtxt(fn, unpack=True, skiprows=1, usecols=0)
        if ray_id_list.size == 0:
            next_ray_id = 0
        elif ray_id_list.size == 1:
            next_ray_id = 1
        else:
            next_ray_id = ray_id_list[-1] + 1
    else:
        outfile = open(fn, 'w')
        outfile.write('ray_id, impact_parameter(kpc), bulk velocity(km/s)[vx, vy, vz], ray_start(kpc)[x, y, z], ray_end(kpc)[x, y, z], sim_center(kpc)[x,y,z] \n')
        outfile.close()
        next_ray_id = 0

    return next_ray_id, fn


def calculate_bulk_los_velocity(bulk_velocity, ray_start_coordinates, ray_end_coordinates):
    ray = np.subtract(ray_end_coordinates, ray_start_coordinates)
    
    bv_mag = np.linalg.norm(bulk_velocity)
    ray_mag = np.linalg.norm(ray)

    # theta is the angle between the ray vector (i.e. line of                                           
    # sight) and the velocity vectors: a dot b = ab cos(theta) 
    cos_theta = np.dot(bulk_velocity, ray) / (bv_mag * ray_mag)
    return bv_mag*cos_theta
    
