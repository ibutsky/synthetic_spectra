from astropy.io import fits
import numpy as np
import glob
import os
import sys
#import pynbody
import yt
import trident
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

def _omass(field, data):
    return data[('Gas', 'Mass')] * data[('Gas', 'OxMassFrac')]

def _femass(field, data):
    return data[('Gas', 'Mass')] * data[('Gas', 'FeMassFrac')]

def _CRPressure(field, data):
    crgamma = 4./3.
    crpressure = (crgamma - 1.) * data[('Gas', 'CREnergy')].d * data[('Gas', 'density')].in_units('g/cm**3').d
    return YTArray(crpressure, 'dyn/cm**2')

def _Pressure(field, data):
    gamma = 5./3.
    m_p = YTQuantity(1.6726219e-24, 'g')  # mass of proton in grams                                                                                                                                                                
    mu = 1.2
    kb = YTQuantity(1.38e-16, 'erg/K') #boltzmann constant cgs                                                                                                                                                                         
    u = data[('Gas', 'Temperature')] * kb / (mu*m_p * (gamma-1.))
    return u * data[('Gas','density')] * (gamma-1.)


def _CRBeta(field, data):
    gamma = 5./3.
    m_p = 1.6726219e-24  # mass of proton in grams                                                                                                                                                                
    mu = 1.2
    kb = 1.38e-16 #boltzmann constant cgs                                                                                                                                                                         
    u = data[('Gas','Temperature')].d * kb / (mu*m_p * (gamma-1.))
    pressure = u * data[('Gas', 'density')].in_units('g/cm**3').d * (gamma-1.)
    #crgamma = 4./3.
 #   crpressure =  (crgamma - 1.) * data[('Gas', 'CREnergy')].d * data[('Density')].in_units('g/cm**3').d                                                                                                         
    return data[('gas', 'cr_pressure')]/ data[('gas', 'pressure')]
    
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


def load_simulation_properties(model, output=3195, ion_list = ['H I', 'O VI', 'Si II', 'Si III', 'Si IV', 'Mg II', 'N V', 'C IV']):
    # note AHF finds this value for P0:  78.13         -239.42          -65.99 
    if model == 'P0':
        ds = yt.load('/Users/irynabutsky/simulations/patient0/pioneer50h243.1536gst1bwK1BH.%06d'%output)
        gcenter = YTArray([-16933.77317667, -12009.28144633,   5305.25448309], 'kpc') # generated with shrink sphere
        bulk_velocity =	YTArray([73.05701672, -239.32976334,  -68.07892736], 'km/s')
        ds.add_field(('gas', 'pressure'), function=_Pressure, sampling_type = 'particle',
                     units = ds.unit_system["pressure"])  
    elif model == 'P0_agncr':
        ds = yt.load('/Users/irynabutsky/simulations/patient0_agncr/pioneer.%06d'%output)
        gcenter = YTArray([-16933.48544591, -12006.24067239,   5307.33807425], 'kpc') # generated with shrink sphere
        bulk_velocity = YTArray([74.98331176, -240.71723683,  -67.77556155], 'km/s')        
        ds.add_field(('gas', 'cr_pressure'), function=_CRPressure, sampling_type = 'particle', 
                     units = ds.unit_system["pressure"])
        ds.add_field(('gas', 'pressure'), function=_Pressure, sampling_type = 'particle',
                     units = ds.unit_system["pressure"])
        ds.add_field(('gas', 'cr_eta'), display_name = ('$P_{\\rm c} / P_{\\rm g}$'),
                     function = _CRBeta, sampling_type = 'particle', units = '')
        
    ds.add_field(('gas', 'O_mass'), function = _omass, sampling_type = 'particle', units = ds.unit_system['mass'])
    ds.add_field(('gas', 'Fe_mass'), function = _femass, sampling_type = 'particle', units = ds.unit_system['mass'])
    trident.add_ion_fields(ds, ions = ion_list)
    
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
    

def calculate_ray_direction(ray_id, ray_data_file = '../../data/random_sightlines.dat'):
    ray_id_list, impact, xi, yi, zi, xf, yf, zf = np.loadtxt(ray_data_file, skiprows = 1, unpack = True)
    
    start_pos = np.array([xi[ray_id], yi[ray_id], zi[ray_id]])
    end_pos   = np.array([xf[ray_id], yf[ray_id], zf[ray_id]])
    return end_pos - start_pos


def load_bulk_velocity(model):
    if model == 'P0':
        return np.array([73.05701672, -239.32976334,  -68.07892736])
    elif model == 'P0_agncr' or model == 'P0agncr':
        return np.array([74.98331176, -240.71723683,  -67.77556155])
