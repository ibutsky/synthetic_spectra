import numpy as np
import os.path
import h5py as h5

import spectrum_generating_tools as spg

# this generates a list of random sightline positions that will be read in by generate_spectra.py
def generate_random_sightline_positions(num_sl = 100, fn = 'random_sightlines.dat', rmin = 10, rmax = 70, ray_len = 500):    
    gcenter = [0, 0, 0] # generating generic coordinates assuming galaxy center is at origin
    ray_id = 0
    
    if os.path.isfile(fn):
        outfile = open(fn, 'a')
        ray_id_list = np.loadtxt(fn, unpack=True, skiprows = 1, usecols=0)
        ray_id = ray_id_list[-1]+1
    else:
        outfile = open(fn, 'w')
        outfile.write('ray_id, impact parameter, ray_start(xyz), ray_end (xyz)\n')

    for i in range(num_sl):
        impact_parameter, ray_start, ray_end = spg.generate_random_ray_coordinates(gcenter, rmin, rmax, ray_len)
        outfile.write("%i %.2f %e %e %e %e %e %e\n"%(ray_id, impact_parameter, ray_start[0], ray_start[1],
                                                     ray_start[2], ray_end[0], ray_end[1], ray_end[2]))
        ray_id += 1

    outfile.close()

rmin = 51
rmax = 300
#rmin = 10
#rmax = 50
generate_random_sightline_positions(rmin = rmin, rmax = rmax)

                


