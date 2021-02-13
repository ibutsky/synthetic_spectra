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
import spectrum_generating_tools as spg

def generate_ray_image_data(field_list, weight_list,
                            model = 'P0', output = 3195, ray_data_file = '',
                            data_loc = '.', ion_list = 'all', redshift = None):

    # load data set with yt, galaxy center, and the bulk velocity of the halo
    ds, gcenter, bulk_velocity = ds, gcenter, bv = spg.load_simulation_properties(model)
    # for annoying reasons... need to convert ray positions to "code units"
    code_unit_conversion = ds.domain_right_edge.d / ds.domain_right_edge.in_units('kpc').d
    print(ray_data_file)
    ray_id_list, impact, bvx, bvy, bvz, xi, yi, zi, xf, yf, zf, cx, cy, cz =\
        np.loadtxt(ray_data_file, skiprows = 1, unpack = True)

    bulk_velocity = YTArray([bvx[0], bvy[0], bvz[0]], 'km/s')
    # set field parameters so that trident knows to subtract off bulk velocity
    ad = ds.all_data()
    ad.set_field_parameter('bulk_velocity', bulk_velocity)
    ad.set_field_parameter('center', gcenter)
    

    width = ds.arr([300., 20., 10], 'kpc')
    for i in [5]:
        # generate the coordinates of the random sightline
        # write ray id, impact parameter, bulk velocity, and start/end coordinates out to file
        h5file = h5.File('%s/ray_image_data_%s_%i_%i.h5'%(data_loc, model, output, ray_id_list[i]), 'w')

        ray_start = ds.arr([xi[i], yi[i], zi[i]], 'kpc')
        ray_end   = ds.arr([xf[i], yf[i], zf[i]], 'kpc')
        ray_direction = ray_end - ray_start
        ray_center = ray_start + 0.5*ray_direction
        normal_vector = ds.arr(gcenter, 'kpc') - ray_center
        normal_vector = np.cross(normal_vector, ray_direction)
        
        for field, weight in zip(field_list, weight_list):
            print(field)
            sys.stdout.flush()
            if weight is not None:
                weight = ('gas', weight)
            if field not in h5file.keys():
                #left_edge = 
                #box = ds.region(ray_center, )
                image = yt.off_axis_projection(ds, center = ray_center, normal_vector = normal_vector,
                width = width, resolution = [1200, 80], item = ('gas', field), weight = weight)
                
                h5file.create_dataset(field, data = image)
                h5file.flush()

        h5file.close()
        print("saved sightline data %i\n"%(i))
            


   
# here's how to actually call this:
ion_list = ['H I', 'O VI', 'Si III']#, 'C II', 'C III', 'C IV', 'Si II', 'Si III', 'Si IV', 'N V']
#ion_list = []
field_list = ['density', 'temperature', 'metallicity', 'velocity_magnitude',
              'velocity_x', 'velocity_y', 'velocity_z', 'O_p5_mass', 'Si_p2_mass', 'H_p0_mass']
           #   'O_p5_number_density', 'Si_p2_number_density', 'H_p0_number_density']
#weight_list = ['density', 'density', 'density', 'density',
#               'density', 'density', 'density', None, None, None]
weight_list = len(field_list)*['density']

data_loc = '../../data/ray_files'
output = 3195

model = 'P0'
ray_data_file = '../../data/unanalyzed_spectra/%s_z0.25_ray_data.dat'%model

generate_ray_image_data(field_list, weight_list, model = model, ion_list = ion_list,
                        ray_data_file = ray_data_file, data_loc = data_loc)


