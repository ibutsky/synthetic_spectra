import matplotlib.pylab as plt
import numpy as np
import sys

import spec_helper_functions as shf

def plot_ion(ax, ion):
    ax.errorbar(wl, flux, yerr=ferr, label = ion)
    ax.axhline(y=0, color='red', linestyle = 'dashed')

    restwave = shf.restwave(ion)
    ax.set_xlim(restwave - 0.5, restwave + 0.5)
    ax.set_ylim(-0.1, 1.5)
    ax.set_xlabel('Wavelength (A)')
    ax.set_ylabel('Normalized Flux')
    ax.legend()


orientation = sys.argv[1]
model = sys.argv[2]
radius = int(sys.argv[3])
if len(sys.argv) > 4:
    ion = sys.argv[4]
else:
    ion = 'all'

fn = shf.spec_folder(orientation, model)+'%ikpc_VPmodel.fits'%(radius)
wl, flux, ferr = shf.load_spec_from_fits(fn)


if ion == 'all':
    ion_list = ['HI', 'CII', 'CIV', 'SiII', 'SiIII', 'SiIV', 'NV']
    nrows = 2
    ncols = int((len(ion_list) + nrows - 1)/ nrows)
    figsize =(3.2*ncols, 3*nrows)
else:
    ion_list = [ion]
    nrows = 1
    ncols = 1
    figsize = (8, 8)


fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize=figsize, sharex = False, sharey = True)

if ion == 'all':
    for i, ion in enumerate(ion_list):
        row = int(i / ncols)
        col = i - row*ncols
        plot_ion(ax[row][col], ion)
else:
    plot_ion(ax, ion)

fig.tight_layout()
plt.show()
