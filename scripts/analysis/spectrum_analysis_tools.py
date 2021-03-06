from astropy.io import fits
import numpy as np
import glob
import os
import sys

def find_spec_impact_parameter(model, redshift, ray_id,  work_dir = '../../data'):
    info_file = '%s/%s_z%0.2f_ray_data.dat'%(work_dir, model, redshift)
    id_list, impact_list = np.loadtxt(info_file, unpack=True, skiprows=1, usecols=(0,1))
    return impact_list[id_list == ray_id]

def extract_spec_info(spec_name, work_dir = '../../data'):
    # note: this assumes a specific naming convention      
    for test_indent in range(3):
        ray_id = spec_name[-1 - test_indent : ]
        if not ray_id.__contains__('_'):
            indent = test_indent

    model = spec_name[8  : -9 - indent]
    print(model)
    redshift = float(spec_name[-6 - indent : -3 - indent])
    ray_id = int(spec_name[-1 - indent :])
    print(ray_id)
    impact = find_spec_impact_parameter(model, redshift, ray_id, work_dir = work_dir)
    return model, redshift, impact, ray_id


def spec_ready_for_analysis(spec_name):
    ready = False
    file1 = '%s/compiledVPoutputs.dat'%(spec_name)
    file2 = '%s/%s_ibnorm.fits'%(spec_name, spec_name)
    if os.path.isfile(file1) and os.path.isfile(file2):
        ready = True
    return ready

def spec_base_filename(model, redshift, ray_id):
    return 'COS-FUV_%s_z%.2f_%i'%(model, redshift, ray_id)

def load_velocity_data(ion, spec, work_dir = '../../data/analyzed_spectra', unanalyzed = False):
    model, redshift, impact, ray_id = extract_spec_info(spec, work_dir = work_dir)
    w0 = restwave(ion, redshift)
    if unanalyzed:
        wl, flux, ferr = load_spec_from_fits('%s/%s.fits'%(work_dir, spec), unanalyzed = True)
    else:
        wl, flux, ferr = load_spec_from_fits('%s/%s/%s_ibnorm.fits'%(work_dir, spec, spec))
    vv = (wl-w0) / w0 * 2.9979e5

    fit = '%s/%s/FitInspection.fits'%(work_dir,spec)
    
    if os.path.isfile(fit):
        wlfit, fluxfit, ferrfit = load_spec_from_fits(fit)
        vvfit = (wlfit-w0) / w0 * 2.9979e5
    else:
        print("ERROR NO VPMODEL FITS FILE")
        wlfit = wl
        fluxfit = flux
        vvfit = vv

    return vv, flux, vvfit, fluxfit, wl, wlfit, w0

def load_spec_from_fits(fn, unanalyzed = False):
    spec = fits.open(fn)
    if unanalyzed:
        data = np.array(spec[1].data[:])
        wl   = np.array([element[0] for element in data])
        flux = np.array([element[1] for element in data])
        err  = np.array([element[2] for element in data])
    else:
        wl = spec['WAVELENGTH'].data
        flux = spec['FLUX'].data
        err = spec['ERROR'].data
    spec.close()
    return wl, flux, err

def all_restwaves(ion, z = 0, infile = 'restwaves_to_use.dat'):
    ion = ion.replace(" ", "")
    ion_names = np.loadtxt(infile, unpack = True, usecols = 0, dtype = 'str')
    restwaves = np.loadtxt(infile, unpack = True, usecols = 1)

    mask = ion_names == ion
    return restwaves[mask]

def restwave(ion, z = 0):
    # note: these were generated by hand                                                                                                                                               
    # using the rest wavelengths generated by                                                                                                                                          
    # the veeper VPmodel fit                                                                                                                                                           
    ion = ion.replace(" ", "")
    if ion == 'HI':
#        restwave = 1215.67
        restwave = 920.963
    elif ion == 'CII':
        restwave = 1334.5323
    elif ion == 'CIII':
        restwave = 977.020
    elif ion == 'CIV':
        restwave = 1548.204
    elif ion == 'NV':
        restwave = 1238.821
    elif ion == 'NIII':
        restwave = 989.799
    elif ion == 'SiII':
        restwave = 1190.4158
    elif ion == 'SiIII':
        restwave = 1206.5
    elif ion == 'SiIV':
        restwave = 1402.7729
    elif ion == 'MgII':
        restwave = 1239.9253
    elif ion == 'OVI':
        restwave = 1031.9261
    else:
        print("%s not recognized as ion"%(ion))
    obs_wavelength = (z+1)*restwave
    return obs_wavelength

def ion_from_restwave(rw):
    if restwave == 1215.67:
        ion  = 'H I'
    elif restwave == 1334.5323:
        ion = 'CII'
    elif restwave == 1548.204:
        ion = 'CIV'
    elif restwave == 1238.821:
        ion = 'NV'
    elif restwave == 1190.4158:
        ion = 'SiII'
    elif restwave == 1206.5:
        ion = 'SiIII'
    elif restwave == 1402.7729:
        ion = 'SiIV'
    else:
        print("%s not a recognized restwave"%(rw))
    return ion

