import yt
from yt.units.dimensions import length
from astropy import constants as const
import sys
import math
import trident
import numpy as np
import os.path

import spec_tools as sp

def make_random_spectrum(model, output, spectrum_directory = '.', \
                 ion_list = 'all', rmin = 10, rmax = 100, ray_len = 500):

    
    ds, gcenter, bulk_velocity = sp.load_simulation_properties(model, output)
    redshift = round(ds.current_redshift, 2)

    kpc_unit = (ds.domain_right_edge[0].d - ds.domain_left_edge[0].d) \
        / (ds.domain_right_edge.in_units('kpc')[0].d - ds.domain_left_edge.in_units('kpc')[0].d)
    ikpc_unit = 1.0 / kpc_unit
    print(ds.domain_right_edge, ds.domain_right_edge.in_units('kpc'), kpc_unit, gcenter)

    
    impact_parameter, ray_start, ray_end = sp.generate_random_ray_coordinates(gcenter, rmin*kpc_unit, rmax*kpc_unit, ray_len*kpc_unit)

    ray_id, ray_fn = sp.get_next_ray_id(model, redshift, spectrum_directory = spectrum_directory)
    ray_outfile = open(ray_fn, 'a')
    ray_outfile.write('%i %.2f %.2f %.2f %.2f %e %e %e %e %e %e %e %e %e\n'%(ray_id, impact_parameter*ikpc_unit, \
                bulk_velocity[0].in_units('km/s').d, bulk_velocity[1].in_units('km/s').d, bulk_velocity[2].in_units('km/s').d,\
                      ray_start[0]*ikpc_unit, ray_start[1]*ikpc_unit, ray_start[2]*ikpc_unit, \
                        ray_end[0]*ikpc_unit,   ray_end[1]*ikpc_unit,   ray_end[2]*ikpc_unit, \
                        gcenter[0]*ikpc_unit, gcenter[1]*ikpc_unit, gcenter[2]*ikpc_unit))
    ray_outfile.close()

    ray = trident.make_simple_ray(ds,
                            start_position = ray_start,
                            end_position = ray_end,
                            lines=ion_list,
                            ftype='gas',
                            # the current redshift of the simulation, calculated above, rounded to two decimal places
                            redshift=redshift,
                            data_filename='%s/ray_%s_%s_%i.h5'%(spectrum_directory, model, output, ray_id))
         
    for instrument in ['COS-G130M', 'COS-G160M']:
        sg = trident.SpectrumGenerator(instrument, line_database = 'prochaska_lines.txt')
        sg.make_spectrum(ray, lines = ion_list)
        sg.apply_lsf()
        sg.add_gaussian_noise(20)
        temp_name = '%s.fits'%(instrument)
        sg.save_spectrum(temp_name, format = 'FITS')
            

    # stitch the spectra together in the format you'll need for veeper
    spec_name = '%s/COS-FUV_%s_z%.2f_%i.fits'%(spectrum_directory, model, redshift, int(ray_id))
    sp.stitch_g130m_g160m_spectra('COS-G130M.fits', 'COS-G160M.fits', spec_name)
    os.remove('COS-G130M.fits')
    os.remove('COS-G160M.fits')
            


         
   
# here's how to actually call this:
ion_list = ['H I', 'O VI', 'C II', 'C III', 'C IV', 'Si II', 'Si III', 'Si IV', 'N V']
sd = '/projects/eot/bafa/data/spectra/tempest'

sd = '../../data/unanalyzed_spectra'
model = 'P0'
output = 3456

for i in range(15):
    make_random_spectrum(model, 3456, ion_list = ion_list, spectrum_directory = sd)

