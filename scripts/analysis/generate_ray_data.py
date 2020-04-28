import yt
from yt import YTArray
from yt.units.dimensions import length
from astropy import constants as const
import sys
import h5py as h5
import math
import trident
import numpy as np
import os.path

sys.path.append('../plotting')
import ion_plot_definitions as ipd

def generate_ray_data(model, output, ray_data_file, data_loc = '.', \
                         ion_list = 'all', redshift = None):
    
    # load data set with yt, galaxy center, and the bulk velocity of the halo
    if model == 'P0':
        ds = yt.load('/nobackup/ibutsky/tmp/pioneer.%06d'%(output))

    trident.add_ion_fields(ds, ions = ion_list)
    # for annoying reasons... need to convert ray positions to "code units"
    code_unit_conversion = ds.domain_right_edge.d / ds.domain_right_edge.in_units('kpc').d
    
    ray_id_list, impact, bvx, bvy, bvz, xi, yi, zi, xf, yf, zf, cx, cy, cz =\
        np.loadtxt('../../data/P0_z0.25_ray_data.dat', skiprows = 1, unpack = True)

    gcenter_kpc = [cx[0], cy[0], cz[0]]  # assuming galaxy center is the same for all sightlines
    gcenter = gcenter_kpc * code_unit_conversion
    bulk_velocity = YTArray([bvx[0], bvy[0], bvz[0]], 'km/s')
    # set field parameters so that trident knows to subtract off bulk velocity
    ad = ds.all_data()
    ad.set_field_parameter('bulk_velocity', bulk_velocity)
    ad.set_field_parameter('center', gcenter)
    

    ray_start_list = np.ndarray(shape=(0, 3))
    ray_end_list = np.ndarray(shape=(0, 3))
    for i in range(len(xi)):
        ray_start_list = np.vstack((ray_start_list, [xi[i], yi[i], zi[i]] * code_unit_conversion))
        ray_end_list   = np.vstack((ray_end_list,   [xf[i], yf[i], zf[i]] * code_unit_conversion))


    # either specify redshift manually, or determine redshift from the redshift of the simulation
    if redshift is None:
        redshift = round(ds.current_redshift, 2)
    print(gcenter, gcenter[0], bulk_velocity)
    
    for i in range(1, 150):
        # generate the coordinates of the random sightline
        # write ray id, impact parameter, bulk velocity, and start/end coordinates out to file
        h5file = h5.File('%s/ray_%s_%i_%i.h5'%(data_loc, model, output, ray_id_list[i]), 'a')
        # generate sightline using Trident
        ray = trident.make_simple_ray(ds,
                            start_position = ray_start_list[i],
                            end_position = ray_end_list[i],
                            lines=ion_list,
                            ftype='gas',
                            field_parameters=ad.field_parameters,
                            # the current redshift of the simulation, calculated above, rounded to two decimal places
                            redshift=redshift)
        ad_ray = ray.all_data()
        # generating the list of all of the 
        field_list = ['y', 'temperature', 'density', 'metallicity', 'dl']
        source_list = [ad, ad, ad, ad, ad_ray]
        unit_list = ['kpc', 'K', 'g/cm**3', 'Zsun', 'cm']
        yt_ion_list = ipd.generate_ion_field_list(ion_list, 'number_density', full_name = False)
        yt_ion_list[0] = 'H_number_density'
        field_list = np.append(field_list, yt_ion_list)
        for j in range(len(yt_ion_list)):
            unit_list.append('cm**-3')
            source_list.append(ad_ray)

        for field,source,unit in zip(field_list, source_list, unit_list):
            if field not in h5file.keys():
                h5file.create_dataset(field, data = source[('gas', field)].in_units(unit))
                h5file.flush()
        h5file.create_dataset('y_lr', data = ad_ray['y'].in_units('kpc'))
        h5file.flush()
        h5file.close()
        print("saved sightline data %i\n"%(i))
            


   
# here's how to actually call this:
ion_list = ['H I', 'O VI', 'C II', 'C III', 'C IV', 'Si II', 'Si III', 'Si IV', 'N V']

data_loc = '../../data/ray_files'
ray_data_file = '../../data/P0_z0.25_ray_data.dat'
model = 'P0'
output = 3195

generate_ray_data(model, output, ray_data_file, ion_list = ion_list, data_loc = data_loc)


