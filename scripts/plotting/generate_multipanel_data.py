#import matplotlib
#matplotlib.use('Agg')

import yt
from yt import YTArray
import trident
import h5py as h5
import numpy as np

import sys
sys.path.append('../analysis')
import spectrum_generating_tools as spg

def _metallicity2(field, data):
    return data[('gas', 'metallicity')].in_units('Zsun')

output = 3195
model = sys.argv[1]
width = 100
view = 'z'


field_list = [('gas', 'density'), ('Gas', 'Temperature'), ('Gas', 'metallicity2'),\
             ('gas', 'O_p5_number_density'), ('gas', 'Si_p2_number_density'),\
              ('gas', 'H_p0_number_density'), ('gas', 'Si_p1_number_density'),\
              ('gas', 'Si_p3_number_density'), ('gas', 'N_p4_number_density'),\
              ('gas', 'C_p3_number_density'), ('gas', 'radial_velocity'),\
              ('gas', 'velocity_relative_z'), ('gas', 'cr_eta')]


weight_field = [('Gas', 'Density'), ('Gas', 'Density'), ('Gas', 'Density'), \
                None, None, None, None, None, None, None, ('Gas', 'Density'), ('Gas', 'Density'), ('Gas', 'Density')]


# load in simulation data and add ion fields
plot_data = h5.File('../../data/simulation_data/multipanel_data_%s_%06d_%i_kpc'%(model, output, width), 'a')

ds, cen, bv = spg.load_simulation_properties(model)

ds.add_field(('Gas', 'metallicity2'), function = _metallicity2, units = 'Zsun', sampling_type = 'particle')
sp = ds.sphere(cen, (500, 'kpc'))
sp.set_field_parameter('bulk_velocity', bv)

left_edge = cen - YTArray([250, 250, 250], 'kpc')
right_edge = cen + YTArray([250, 250, 250], 'kpc')

box = ds.region(cen, left_edge, right_edge)
box.set_field_parameter('bulk_velocity', bv)


# set up projection plots for fields that are weighted and unweighted
#del plot_data['radial_velocity']
for i in range(len(field_list)):
    print(field_list[0])
    field = field_list[i][1]
    dset = '%s_%s'%(field,view)
    if dset not in plot_data.keys():
        proj = yt.ProjectionPlot(ds, view, field, weight_field = weight_field[i], width=(width, 'kpc'), center = cen, data_source = box)
        proj_frb =  proj.data_source.to_frb((width, 'kpc'), 800)
    
        if field_list[0].__contains__('velocity'):
            plot_data.create_dataset(dset, data = np.array(proj_frb[field_list[i]].in_units('km/s')))
        if field_list[0] == 'cr_eta' and model == 'P0':
            continue
        else:
            plot_data.create_dataset(dset, data = np.array(proj_frb[field_list[i]]))
        plot_data.flush()



