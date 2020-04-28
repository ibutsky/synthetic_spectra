import os
import sys
import numpy as np
import glob

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import ion_plot_definitions as ipd

import plotting_tools as pt 

sys.path.append('../analysis')
import spectrum_analysis_tools as spa

import seaborn as sns
sns.set_style("ticks",{'axes.grid': True, "ytick.major.size": 0.1,
                "ytick.minor.size": 0.05,
                'grid.linestyle': '--'
            })

ion_list = ['H I', 'Si II', 'C III', 'O VI']
color_list = ['gray', 'black', 'black', 'red']

spec_list = glob.glob('../../data/analyzed_spectra/COS-FUV*')
for spec_path in spec_list:
    plot_name = '%s/velocity.png'%(spec_path)
    if not os.path.isfile(plot_name):
        spec = os.path.basename(spec_path)
        print(spec)
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 3.8))
        for ion, color in zip(ion_list, color_list):
            vv, flux, vvfit, fluxfit,  wl, wlfit, w0 = spa.load_velocity_data(ion, spec)
            #ax.plot(vv, flux, color = 'black', alpha = 0.2)
            if os.path.isfile('%s/FitInspection.fits'%(spec_path)):
                ax.plot(vvfit, fluxfit, color = 'black', linewidth = 4, alpha = 0.7)
                ax.plot(vvfit, fluxfit, color = color, label = ion, linewidth = 2, alpha = 0.8)
        ax.set_xlim(-300, 300)
        ax.set_ylim(-0.05, 1.3)
        ax.set_xlabel('Relative Velocity (km/s)')
        ax.set_ylabel('Normalized Flux')
        ax.legend(frameon = True)
        fig.tight_layout()
        plt.savefig(plot_name, dpi = 300)

