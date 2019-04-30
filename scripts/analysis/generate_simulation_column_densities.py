import yt
import trident
import h5py as h5
import numpy as np

import ion_plot_definitions as ipd

def generate_sim_column_densities(sim, ion_list, res = 800):
    cdens_file = h5.File('../data/simulated_ion_column_densities_%s.h5'%(sim), 'a')
    
    ds = yt.load('/Users/irynabutsky/Work/galaxy/g160_torB_%s/DD2600'%(sim))
    trident.add_ion_fields(ds, ion_list)

    v, center = ds.h.find_max(("gas", "density")) 

    
    field_list = ipd.generate_ion_field_list(ion_list, 'number_density')    
    print(field_list)
  
    width = yt.YTQuantity(200, 'kpc')
    px, py = np.mgrid[-width/2:width/2:res*1j, -width/2:width/2:res*1j]
    radius = (px**2.0 + py**2.0)**0.5
    
    if "radius" not in cdens_file.keys():
        cdens_file.create_dataset("radius", data = radius.ravel())
        
    for axis in ['x', 'y', 'z']:
        frb = ipd.make_projection(ds, axis, field_list, center, width)
        
        for i, ion in enumerate(ion_list):
            dset = "%s_%s" % (ion.replace(" ", ""), axis)
            if dset not in cdens_file.keys():
                cdens_file.create_dataset(dset, data=frb[field_list[i]].ravel())
                cdens_file.flush()



ion_list = ['H I', 'C II', 'C IV', 'Si II', 'N V', 'Si III', 'Si IV', 'Mg II']

generate_sim_column_densities('ad', ion_list)
generate_sim_column_densities('s4', ion_list)

