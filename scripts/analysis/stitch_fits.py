from astropy.io import fits
import numpy as np
import glob
import sys
import os


def stitch_g130m_g160m_spectra(fn_g130, fn_g160, fn_combined):
    temp130 = fits.open(fn_g130)
    temp160 = fits.open(fn_g160)
    data = []
    for item in temp130[1].columns:
        if item.name in ['wavelength', 'flux', 'flux_error']:
            # this iterates through each of the arrays in the extension by grabbing names from the header           
            data1     = temp130[1].data[item.name] # get data from file 1        
            data2     = temp160[1].data[item.name]            
            new_array = np.concatenate([data1, data2]) # make the new array composed of both parts    
            data.append(new_array)

    c1 = fits.Column(name = 'WAVELENGTH', array = data[0], format = 'E')
    c2 = fits.Column(name = 'FLUX', array = data[1], format = 'E')
    c3 = fits.Column(name = 'FLUXERR', array = data[2], format = 'E')
    fits_combined = fits.BinTableHDU.from_columns([c1, c2, c3])
    fits_combined.writeto(fn_combined)

def stitch_individual_spectrum(view, model, time, impact):
# Function: Stitches together G130M and G160M spectra generated with Trident into a single spectrum 
# Parameters: all parameters are used to specify the name of the spectra files.  
#             view   - Orientation from which the spectrum was generated. Possible options:
#                       'face', 'edge_theta0', 'edge_theta1.0', 'edge_theta1.5'
#             model  - The cosmic ray transport approximation in the simulation. Possible options: 
#                        'anisd' (anisotropic diffusion) or 'stream' (anisotropic streaming)
#             time   - The galaxy's age at the time the spectrum was generated. 
#             impact - The impact parameter from the center of the galaxy.

    base_dir = '../../data/unanalyzed_spectra/'
    base     = '_%s_%s_%sGyr_r%skpc.fits'%(view, model, time, r)      
    fn_g130     = base_dir + 'COS-G130M'+ base
    fn_g160     = base_dir + 'COS-G160M'+ base
    fn_combined = base_dir + 'COS-FUV'          + base

    stitch_g130m_g160m_spectra(fn_g130, fn_g160, fn_combined)
    # move original files to trash folder; this will be cleaned out periodically
    trash_dir = '../../data/trash/'
    os.rename(fn_g130, trash_dir + 'COS-G130M' + base)
    os.rename(fn_g160, trash_dir + 'COS-G160M' + base)

def stitch_all_spectra(dir = '../../data/unanalyzed_spectra/'):
# Function: stitch all of the G130M and G160M spectra files in a given directory
#

    os.chdir(dir)
    g130_files = glob.glob('COS-G130M*')
#    g160_files = glob.glob('COS-G160M*')

    # note: there should be an equal number of g130 and g160 files
    # either way, this only does anything if both the g130 spectrum and 
    # the g160 spectrum are present. 
    # Arbitrarily choosing to loop through g130 files
    for g130_file in g130_files:
        base = g130_file[9:]
        g160_file = 'COS-G160M'+base
        fn_combined = 'COS-FUV'+base
        if os.path.isfile(g160_file) and not os.path.isfile(fn_combined):
            stitch_g130m_g160m_spectra(g130_file, g160_file, fn_combined)




#view = sys.argv[1]
#model = sys.argv[2]
#time = sys.argv[3]
#r = sys.argv[4]

#stitch_individual_spectrum(view, model, time, r)

stitch_all_spectra()
