import matplotlib.pylab as plt
import numpy as np
import sys

import spec_helper_functions as shf
import eqwrange as eqw


#orientation = sys.argv[1]
#model = sys.argv[2]
#radius = int(sys.argv[3])
#if len(sys.argv) > 4:
#    ion = sys.argv[4]
#else:
#    ion = 'all'

orientation = 'face'
model = 'stream'
radius = 20
ion = 'CII'

fn = shf.spec_folder(orientation, model)+'%ikpc_VPmodel.fits'%(radius)
wl, flux, ferr = shf.load_spec_from_fits(fn)    

ion_names = np.loadtxt('../dla.lst', unpack=True, skiprows = 1, usecols=(1), dtype = 'str')
ion_wls, ion_fs = np.loadtxt('../dla.lst', unpack=True, skiprows = 1,usecols=(0, 3))
ion_wls_int = ion_wls.astype(int)

w0 = shf.restwave(ion) # REST wavelength for rest-frame EW; observed wavelength for observed EW
ion_mask  = (ion_names == ion) 
wl_mask =  (ion_wls_int[ion_mask] == int(w0))
f0 = ion_fs[ion_mask][wl_mask][0]

print('*****')
print(w0, f0)
vrange = (-100,100)


range_mask = ((wl > (w0 - 1)) & (wl < (w0 + 1)))
EW, colm, flag, velcent, velwidth = eqw.eqwrange(wl, flux, ferr, vrange, w0, f0)
print(EW, colm, flag, velcent, velwidth)
