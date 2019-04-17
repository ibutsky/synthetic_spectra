import yt
from yt.units.dimensions import length 

import sys
import math
import trident
import numpy as np
import os.path


def make_spectra(model_name, output, orientation, theta_list = [0], spectrum_directory = '.', \
                 ion_list = ['H I', 'O VI'], rmin = 10, rmax = 100, delta_r = 10, ray_length = 1000):
    # model_name:  The name of the GM model... i.e. "P0", "GM1", "GM2"..
    #              this convention is really up to you. These model name
    #              tags will appear throughout the analysis process.
    # orientation: The direction of the line-of-sight of the ray. Options: 'x', 'y', 'z'. 
    #              Default: 'x'. 
    #              Note: this just assumes that 'z' is face-on.
    # theta_list:  A list of angles in which to generate light-rays (in radians). 
    #              Default: theta_list = [0]. Note: theta = 0 lies along the horizontal axis, to the right.
    # spectrum_directory: The directory where the generated spectra will be
    #                     saved. Default: the current directory.
    # ion_list:   The names of the ions to include in your spectra. Default:
    #             H I and O VI
    # rmin:       The minimum impact parameter (in kpc). Default: 10  kpc
    # rmax:       The maximum impact parameter (in kpc). Default: 100 kpc
    # delta_r:    The step-size (in kpc) between rmin and rmax; Default: 10 kpc
    # ray_length: The total length of the ray, in kpc.  If the ray is too long, and moves into
    #             a region with no gas particles, you'll get an error.

    # WARNING: you probably need to change the definition of 'fn' to match your directory structure
    if model_name == 'P0':
        sim_folder = 'pioneer50h243.1536g1bwK1BH'
        sim_name = 'pioneer50h243.1536gst1bwK1BH'
    else:
        sim_name = 'pioneer50h243GM1.1536gst1bwK1BH'%(model_name)
        sim_name = sim_folder

    fn = '/nobackupp2/nnsanche/%s/%s.%06d'%(sim_folder, sim_name, output)
    print(fn)
    fn = '/nobackup/ibutsky/tmp/pioneer.003456'
    ds = yt.load(fn)
	
    # this is how you would add the ions as "fields", but you don't need this
    # to generate spectra 
    # trident.add_ion_fields(ds, ions= ion_list)

    # array of impact parameters
    impact_kpc = np.arange(rmin, rmax+delta_r, delta_r)
	
    time = ds.current_time.in_units('Gyr')
    redshift = ds.current_redshift

    #  find the galaxy center
    v, cen = ds.h.find_max(("gas", "density"))
	
    # get rid of the units, only keep the number value
    # note gcenter is in 'code_length' units, not kpc
    gcenter = cen.d
    print('domain_center', ds.domain_center)
    print(gcenter, cen.in_units('kpc'))
    sys.stdout.flush()
	
    # definte the list of angles if one hasn't been defined
    if theta_list is None:
        if orientation == 'face':
            theta_list = [0]
        elif orientation == 'edge':
            theta_list = [0, np.pi/3., np.pi/2.]
	
    for theta in theta_list:
        for r_kpc in impact_kpc:
            r = (r_kpc / ds.length_unit).d
            # calculates the position of the ray's start and end points
            # depends on the orientation
            # goes through the simulation box along the line of sight
            
            half_ray = (ray_length / 2.0 / ds.length_unit).d
            # first, we define this assuming we're looking into the x direction
            ray_start = [-half_ray, r*np.cos(theta), r*np.sin(theta)]
            ray_end   = [half_ray,  r*np.cos(theta), r*np.sin(theta)]

            # if the orientation is in the y or z directions, we shift the indices
            if orientation == 'y':
                ray_start = np.roll(ray_start, 1)
                ray_end   = np.roll(ray_end, 1)
            elif orientation == 'z':
                ray_start = np.roll(ray_start, 2)
                ray_end   = np.roll(ray_end, 2)

            ray_start += gcenter
            ray_end   += gcenter

            print(ray_start, ray_end)

            # this makes a 'light_ray' object and adds the absorption lines of
            # the elements in ion_list
            # Outputs some information in a ray___.h5 file. You probably
            # won't need this file, but if you parallelize this and have
            # multiple programs try to open the same ray.h5 file, you'll run
            # into problems, hence they're all named uniquely.
            ray = trident.make_simple_ray(ds,
                            start_position = ray_start,
                            end_position = ray_end,
                            lines=ion_list,
                            ftype='gas',
                            # the current redshift of the simulation, calculated above
                            redshift=redshift,
                            data_filename='ray_%s_%i.h5'%(model_name, int(redshift)))
         
            # this actually makes the absorption spectrum for both the G130M and
            # G160 M instruments and saves them as two separate fits files
										  
            for instrument in ['COS-G130M', 'COS-G160M']:
                # generates the spectrum using the specified COS instrument
               sg = trident.SpectrumGenerator(instrument)
               # adds absorption lines
               sg.make_spectrum(ray, lines = ion_list)
               # adds instrument-specific line spread-function
               sg.apply_lsf()
               # adds noise with a signal-to-noise ratio of 10
               sg.add_gaussian_noise(10)

               spec_name = '%s_%s_%s_theta%.1f_z%.2f_r%ikpc.fits'%(instrument, model_name, orientation, theta, redshift, int(r_kpc))

               sg.save_spectrum(spec_name, format = 'FITS')

         
   
# here's how to actually call this:
make_spectra('P0', 3456, 'x')


