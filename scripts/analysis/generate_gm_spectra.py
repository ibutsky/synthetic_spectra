import yt
from yt.units.dimensions import length 

import sys
import math
import trident
import numpy as np
import os.path


def make_spectra(model_name, output, orientation, theta_list = None, \
				spectrum_directory = '.', ion_list = ['H I', 'O VI']):
	# model_name:  The name of the GM model... i.e. "P0", "GM1", "GM2"..
	#              this convention is really up to you. These model name
	#              tags will appear throughout the analysis process.
	# orientation: The orientation of the galaxy. Options are: 'face' and 'edge'.
	#              Note: this just assumes that 'z' is face-on.
	# theta_list:  A list of custom angles in which to generate light-rays.
	#              Default: None (the program makes the list).
	# spectrum_directory: The directory where the generated spectra will be
	#                     saved. Default: the current directory.
	# ion_list:   The names of the ions to include in your spectra. Default:
	#             H I and O VI
	
	# change
    fn = ''
    ds = yt.load(fn)
	
	# this is how you would add the ions as "fields", but you don't need this
	# to generate spectra
	#    trident.add_ion_fields(ds, ions= ion_list)
   
   # array of impact parameters
    impact_kpc = np.linspace(rmin, rmax, delta_r)
	
	time = ds.current_time.in_units('Gyr')
	redshift = ds.current_redshift

	#  find the galaxy center
    v, cen = ds.h.find_max(("gas", "density"))
	
	# get rid of the units, only keep the number value
	# note gcenter is in 'code_length' units, not kpc
    gcenter = [cen[1].d, cen[2].d, cen[3].d]
	
	
	# definte the list of angles if one hasn't been defined
	if theta_list is None:
		if orientation == 'face':
			theta_list = [0]
		elif orientation == 'edge':
			theta_list = [0, np.pi/3., np.pi/2.]
	
    for theta in theta_list:
        for i in range(len(impact)):
            r_kpc = impact_kpc[i]
			r = (r_kpc / ds.length_unit).d
			# calculates the position of the ray's start and end points
			# depends on the orientation
			# goes through the entire simulation box along the line of sight
			#  (0 to 1 in code units)
            if orientation == 'edge':
                ray_start = [0, gcenter[1]+r*np.sin(theta), gcenter[2]+r*np.cos(theta)]
                ray_end = [1, gcenter[1]+r*np.sin(theta), gcenter[2]+r*np.cos(theta)]

            elif orientation == 'face':
                ray_start = [gcenter[1]+r*np.sin(theta), gcenter[2]+r*np.cos(theta), 0]
                ray_end = [gcenter[1]+r*np.sin(theta), gcenter[2]+r*np.cos(theta), 1]

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
							# the current redshift of the simulation
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
               if orientation == 'edge':
				   spec_name = '../data/spectra/rotation_test/%s_edge_theta%.1f_\
				   	%s_%.2fGyr_r%ikpc.fits'%(instrument, model_name, theta, redshift, r_kpc))
               elif orientation == 'face':
                   spec_name = '../data/spectra/rotation_test/%s_face_%s_z%.2f\
				   _r%ikpc.fits'%(instrument, model_name, redshift, r_kpc)
               sg.save_spectrum(spec_name, format = 'FITS')

         
   
# here's how to actually call this:
make_spectra('P0', 2240,
for galaxy in ['g160_torB_s4', 'g160_torB_ad']:
    make_spectra(galaxy, 2240, 0.2, 'face', [0])
    make_spectra(galaxy, 2240, 0.2, 'edge', [0, np.pi/3., np.pi/2.])
   # make_spectra(galaxy, 2600, 'edge', [0, np.pi/4., np.pi/3., np.pi/2.])


