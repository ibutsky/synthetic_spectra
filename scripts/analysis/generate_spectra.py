import yt
from yt.units.dimensions import length 

import sys
import math
import trident
import numpy as np
import os.path

def make_spectra(galaxy, output, redshift, orientation, theta_list):
    
    loc = '/scratch/eot/butsky/'
    ds = yt.load(loc+galaxy+'/DD%04d/DD%04d'%(int(output),int(output)))

#    trident.add_ion_fields(ds, ions=['H I', 'Mg II', 'Si III', 'O VI'])
   
    if galaxy == 'g160_torB_s4':
        crname = 'stream'
    elif galaxy == 'g160_torB_ad':
        crname = 'anisd'

 #   line_list =['Si II', 'Mg II', 'C II', 'N II', 'Si III', 'C III', 'N III', \
  #                                 'C IV', 'Si IV', 'N V', 'O VI', 'H I', 'H II']
        
    ikpc_unit = 0.00076294
#    impact = np.linspace(0.008, 0.15, 100) # desired impact parameter in code_length  
    impact_kpc = np.linspace(10, 100, 10)
    impact = impact_kpc*ikpc_unit
    rmin = 0.008
    rmax = 0.15
    time = int(output)/200. 

#    rmin = yt.YTQuantity(10, 'kpc')
#    rmax = yt.YTQuantity(200, 'kpc')


    v, cen = ds.h.find_max(("gas", "density")) # find the center to keep the galaxy at the center of all the images;
    sp = ds.sphere(cen, (30.0, "kpc"))
    cen2 = sp.quantities.center_of_mass(use_gas=True, use_particles=False).in_units("kpc")
    sp2 = ds.sphere(cen2, (1.0, "kpc"))
    cen3 = sp2.quantities.max_location(("gas", "density"))
    gcenter = [cen3[1].d, cen3[2].d, cen3[3].d]


    for theta in theta_list:
        if orientation == 'edge':
            column_outfile = open('../data/spectra/rotation_test/column_edge_theta'+str(theta)[:3]+'_'+crname+'_'+str(time)+'Gyr_r', 'w')
        elif orientation == 'face':
            column_outfile = open('../data/spectra/rotation_test/column_face__'+crname+'_'+str(time)+'Gyr_r', 'w')

            column_outfile.write("# radius(kpc)  H I, Mg II, SiIII, O VI (column densities in units of 1/cm^2\n")

        for i in range(len(impact)):
            r = impact[i]
            r_kpc = impact_kpc[i]
            
            if orientation == 'edge':
                ray_start = [0, gcenter[1]+r*np.sin(theta), gcenter[2]+r*np.cos(theta)]
                ray_end = [1, gcenter[1]+r*np.sin(theta), gcenter[2]+r*np.cos(theta)]

            elif orientation == 'face':
                ray_start = [gcenter[1]+r*np.sin(theta), gcenter[2]+r*np.cos(theta), 0]
                ray_end = [gcenter[1]+r*np.sin(theta), gcenter[2]+r*np.cos(theta), 1]

            ray = trident.make_simple_ray(ds,
                            start_position = ray_start,
                            end_position = ray_end,
                            lines='all',
                            ftype='gas',
                            redshift=redshift,
                            data_filename='ray_%s_%i.h5'%(galaxy, int(output)))
         
            for instrument in ['COS-G130M', 'COS-G160M']:
               sg = trident.SpectrumGenerator(instrument)
               sg.make_spectrum(ray, lines= 'all')
               sg.apply_lsf()
               sg.add_gaussian_noise(10)
               time = int(output)/200.
               if orientation == 'edge':
                   spec_name = '../data/spectra/rotation_test/'+instrument+'_edge_theta'+str(theta)[:3]+'_'+crname+'_'+str(time)+'Gyr_r'+str(int(r_kpc))+'kpc.fits'
               elif orientation == 'face':
                   spec_name = '../data/spectra/rotation_test/'+instrument+'_face_'+crname+'_'+str(time)+'Gyr_r'+str(int(r_kpc))+'kpc.fits'
               sg.save_spectrum(spec_name, format = 'FITS')

         
   
#for galaxy in ['g160_torB_s4', 'g160_torB_ad']:
for galaxy in ['g160_torB_s4', 'g160_torB_ad']:
    make_spectra(galaxy, 2240, 0.2, 'face', [0])
    make_spectra(galaxy, 2240, 0.2, 'edge', [0, np.pi/3., np.pi/2.])
   # make_spectra(galaxy, 2600, 'edge', [0, np.pi/4., np.pi/3., np.pi/2.])


