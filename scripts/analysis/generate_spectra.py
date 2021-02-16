import yt
from yt import YTArray
import trident
import numpy as np
import os.path
import sys

import spectrum_generating_tools as spg

def make_random_spectrum(model, output = 3195, spectrum_directory = '../../data/unanalyzed_spectra', \
                         ion_list = 'all', redshift = None):

    ds, gcenter, bulk_velocity = spg.load_simulation_properties(model)
    
    ad = ds.all_data()
    ad.set_field_parameter('bulk_velocity', bulk_velocity)
    ad.set_field_parameter('center', gcenter)

    if redshift is None:
        redshift = round(ds.current_redshift, 2)
    print(gcenter, gcenter[0], bulk_velocity)
        
    ray_id_list, impact_list, xs, ys, zs, xe, ye, ze = np.loadtxt('../../data/random_sightlines.dat', skiprows=1, unpack = True)    
    spectrum_data_filename = '%s/%s_z%0.2f_ray_data.dat'%(spectrum_directory, model, redshift)
    total_column_filename = '%s/%s_total_column.dat'%(spectrum_directory, model)
    if not os.path.isfile(spectrum_data_filename):
        spectrum_data_outfile = open(spectrum_data_filename, 'w')
        spectrum_data_outfile.write('ray_id, impact_parameter(kpc), bulk velocity(km/s)[vx, vy, vz], ray_start(kpc)[x, y, z], \
        ray_end(kpc)[x, y, z], sim_center(kpc)[x,y,z] \n')
    else:
        spectrum_data_outfile = open(spectrum_data_filename, 'a')
    if not os.path.isfile(total_column_filename):
        col_outfile = open(total_column_filename, 'w')
        col_outfile.write('#ray_id, impact_parameter(kpc), OVI col, SiIII col, HI col\n')
    else:
        col_outfile = open(total_column_filename, 'a')
        
    for i, ray_id in enumerate(ray_id_list):
        spec_name = '%s/COS-FUV_%s_z%.2f_%04d.fits'%(spectrum_directory, model, redshift, int(ray_id))
        if not os.path.isfile(spec_name):
            ray_start = ds.arr(gcenter.in_units('kpc').d + np.array([xs[i], ys[i], zs[i]]), 'kpc')
            ray_end   =	ds.arr(gcenter.in_units('kpc').d + np.array([xe[i], ye[i], ze[i]]), 'kpc')
            spectrum_data_outfile.write('%i %.2f %.2f %.2f %.2f %e %e %e %e %e %e %e %e %e\n'%(ray_id, impact_list[i],
                                         bulk_velocity[0], bulk_velocity[1], bulk_velocity[2], ray_start[0], 
                                         ray_start[1], ray_start[2], ray_end[0], ray_end[1], ray_end[2],
                                         gcenter[0].astype(float), gcenter[1].astype(float), gcenter[2].astype(float)))
            spectrum_data_outfile.flush()

            ray = trident.make_simple_ray(ds,
                            start_position = ray_start,
                            end_position = ray_end, 
                            lines=ion_list,
                            ftype='gas',
                            fields           = ['temperature', 'density', 'metallicity', 'velocity_x', 'velocity_y', 'velocity_z'],
                            field_parameters=ad.field_parameters,
                            # the current redshift of the simulation, calculated above, rounded to two decimal places
                            redshift=redshift,
                            data_filename='%s/ray_%s_%s_%i.h5'%(spectrum_directory, model, output, ray_id))

            rd = ray.all_data()
            dl = rd['dl']
            ocol = np.cumsum(rd['O_p5_number_density']*dl)[-1]
            sicol = np.cumsum(rd['Si_p2_number_density']*dl)[-1]
            hcol = np.cumsum(rd['H_p0_number_density']*dl)[-1]

            col_outfile.write("%i %.2f %.2f %.2f %.2f\n"%(ray_id, impact_list[i], np.log10(ocol), np.log10(sicol), np.log10(hcol)))
            col_outfile.flush()
            
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
    col_outfile.close()

   
# here's how to actually call this:
ion_list = ['H I', 'O VI', 'C II', 'C III', 'C IV', 'Si II', 'Si III', 'Si IV', 'N V', 'Mg II']

sd = '../../data/unanalyzed_spectra'

#model = 'P0'
model = sys.argv[1]

for i in range(15):
    make_random_spectrum(model, ion_list = ion_list, spectrum_directory = sd)


