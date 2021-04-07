import yt
from yt import YTArray
import trident
import h5py as h5
import numpy as np

import sys
sys.path.append('../plotting')

import ion_plot_definitions as ipd
import spectrum_generating_tools as spg

def generate_sim_column_densities(sim, ion_list, res = 800):
    #cdens_file = h5.File('../../data/simulated_ion_column_densities_%s.h5'%(sim), 'w')
    cdens_file = h5.File('test.h5', 'w')

    ds, center, bv = spg.load_simulation_properties(sim)
    print(center)
    center = ds.arr(center, 'kpc')
        
    trident.add_ion_fields(ds, ion_list)
    
    field_list = ipd.generate_ion_field_list(ion_list, 'number_density', model = 'P0')

    width = yt.YTQuantity(110, 'kpc')
    px, py = np.mgrid[-width/2:width/2:res*1j, -width/2:width/2:res*1j]
    radius = (px**2.0 + py**2.0)**0.5
    
    if "radius" not in cdens_file.keys():
        cdens_file.create_dataset("radius", data = radius.ravel())

    # note: something weird going on in x-axis with P0
    for axis in ['x', 'y', 'z']:
        frb = ipd.make_projection(ds, axis, field_list, center, width)
        
        for i, ion in enumerate(ion_list):
            dset = "%s_%s" % (ion.replace(" ", ""), axis)
            if dset not in cdens_file.keys():
                cdens_file.create_dataset(dset, data=frb[field_list[i]].ravel())
                cdens_file.flush()



ion_list = ['O VI','Si III']#, 'Si II', 'Si IV', 'C III']#, 'C II', 'C IV', 'Si II', 'N V', 'Si III', 'Si IV', 'Mg II']

#generate_sim_column_densities('ad', ion_list)
#generate_sim_column_densities('s4', ion_list)
for model in ['P0']:#, 'P0_agncr']:
    generate_sim_column_densities(model, ion_list)

