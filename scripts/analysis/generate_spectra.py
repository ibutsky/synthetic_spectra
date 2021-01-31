import yt
from yt import YTArray
import trident
import numpy as np
import os.path
import sys

import spectrum_generating_tools as spg

def make_random_spectrum(model, output = 3195, spectrum_directory = '../../data/unanalyzed_spectra', \
                         ion_list = 'all', redshift = None):

    if model == 'P0':
        ds = yt.load('/Users/irynabutsky/simulations/patient0/pioneer.%06d'%output)
        gcenter = YTArray([-1.693207e4, -1.201068e4, 5.303337e3], 'kpc')
        bulk_velocity = YTArray([72.78, -248.83, -63.49], 'km/s')
    if model == 'P0_agncr':
        ds = yt.load('/Users/irynabutsky/simulations/patient0_agncr/pioneer.%06d'%output)
        gcenter = YTArray([-1.695524e4, -1.199842e4, 5.309064e3], 'kpc')
        bulk_velocity = YTArray([74.00, -249.27, -63.50], 'km/s')

    ad = ds.all_data()
    ad.set_field_parameter('bulk_velocity', bulk_velocity)
    ad.set_field_parameter('center', gcenter)

    if redshift is None:
        redshift = round(ds.current_redshift, 2)
    print(gcenter, gcenter[0], bulk_velocity)
    gcenter = np.array(gcenter.d)
    
    kpc_unit = (ds.domain_right_edge[0].d - ds.domain_left_edge[0].d) \
        / (ds.domain_right_edge.in_units('kpc')[0].d - ds.domain_left_edge.in_units('kpc')[0].d)
    ikpc_unit = 1.0 / kpc_unit
    print(kpc_unit)
    
    ray_id_list, impact_list, xs, ys, zs, xe, ye, ze = np.loadtxt('../../data/random_sightlines.dat', skiprows=1, unpack = True)    
    spectrum_data_filename = '%s/%s_z%0.2f_ray_data.dat'%(spectrum_directory, model, redshift)
    if not os.path.isfile(spectrum_data_filename):
        spectrum_data_outfile = open(spectrum_data_filename, 'w')
        spectrum_data_outfile.write('ray_id, impact_parameter(kpc), bulk velocity(km/s)[vx, vy, vz], ray_start(kpc)[x, y, z], \
        ray_end(kpc)[x, y, z], sim_center(kpc)[x,y,z] \n')
    else:
        spectrum_data_outfile = open(spectrum_data_filename, 'a')
        
    for i, ray_id in enumerate(ray_id_list):
        spec_name = '%s/COS-FUV_%s_z%.2f_%04d.fits'%(spectrum_directory, model, redshift, int(ray_id))
        if not os.path.isfile(spec_name):
            ray_start = YTArray(gcenter + np.array([xs[i], ys[i], zs[i]]), 'kpc')
            ray_end   =	YTArray(gcenter + np.array([xe[i], ye[i], ze[i]]), 'kpc')
            spectrum_data_outfile.write('%i %.2f %.2f %.2f %.2f %e %e %e %e %e %e %e %e %e\n'%(ray_id, impact_list[i],
                                         bulk_velocity[0], bulk_velocity[1], bulk_velocity[2], ray_start[0], 
                                         ray_start[1], ray_start[2], ray_end[0], ray_end[1], ray_end[2],
                                         gcenter[0].astype(float), gcenter[1].astype(float), gcenter[2].astype(float)))
            spectrum_data_outfile.flush()

            ray = trident.make_simple_ray(ds,
                            start_position = ds.arr(ray_start.d*kpc_unit, 'unitary'),
                            end_position = ds.arr(ray_end.d*kpc_unit, 'unitary'),
                            lines=ion_list,
                            ftype='gas',
                            field_parameters=ad.field_parameters,
                            # the current redshift of the simulation, calculated above, rounded to two decimal places
                            redshift=redshift,
                            data_filename='%s/ray_%s_%s_%i.h5'%(spectrum_directory, model, output, ray_id))
         
            for instrument in ['COS-G130M', 'COS-G160M']:
                sg = trident.SpectrumGenerator(instrument, line_database = 'pyigm_line_list.txt')
                sg.make_spectrum(ray, lines = ion_list)
                sg.apply_lsf()
                sg.add_gaussian_noise(10)
                temp_name = '%s_%s.fits'%(instrument, model)
                sg.save_spectrum(temp_name, format = 'FITS')
            

            # stitch the spectra together in the format you'll need for veeper
            spg.stitch_g130m_g160m_spectra('COS-G130M_%s.fits'%model, 'COS-G160M_%s.fits'%model, spec_name)
            os.remove('COS-G130M_%s.fits'%model)
            os.remove('COS-G160M_%s.fits'%model)
            
    spectrum_data_outfile.close()

   
# here's how to actually call this:
ion_list = ['H I', 'O VI', 'C II', 'C III', 'C IV', 'Si II', 'Si III', 'Si IV', 'N V', 'Mg II']

sd = '../../data/unanalyzed_spectra'

#model = 'P0'
model = sys.argv[1]

for i in range(15):
    make_random_spectrum(model, ion_list = ion_list, spectrum_directory = sd)


