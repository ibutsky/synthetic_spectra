import matplotlib.pyplot as plt

import numpy as np
import h5py as h5
import sys
from matplotlib.colors import LogNorm

import seaborn as sns
sns.set_style("white",{'font.family':'serif', 'axes.grid': True, "ytick.major.size": 0.1,
                "ytick.minor.size": 0.05,
                'grid.linestyle': '--'})

import ion_plot_definitions as ipd
    
def plot_multipanel(ion_list, sim):
    if sim == 'ad':
        model = 'anisd'
    elif sim == 's4':
        model = 'stream'

    nrows = 2
    ncols = int((len(ion_list)+1)/2)
    fig, ax = plt.subplots(nrows = nrows, ncols =ncols , figsize = (4.2*ncols, 4*nrows),sharex = True, sharey = False)

    fn = '../data/simulated_ion_column_densities_%s.h5'%(sim)
    
    for i, ion in enumerate(ion_list):
        row = int(i/ncols)
        col = int(i - ncols*row)

        ylims = ipd.return_ylims(ion)
        rmax = 100

        r_arr, cdens_arr = ipd.load_r_cdens(fn, ion)
        im =  ipd.plot_hist2d(ax[row][col], r_arr, cdens_arr, rmax,  ylims, vmin = 1e-3)
            
        xscatter, yscatter, yerr, flagscatter = ipd.load_veeper_ion_column(ion.replace(" ", ""), model)
#        print(impact, col)
        yscatter = 10**yscatter
        yerr = 10**yerr
#        ax.set_yscale('log')
        
        sat = flagscatter > 8
        ax[row][col].scatter(xscatter[sat], yscatter[sat],  marker = "^", s = 40, edgecolor = 'black', facecolor = 'orange')
        unsat = flagscatter <= 8
        ax[row][col].scatter(xscatter[unsat], yscatter[unsat],  marker = "s", s = 40, edgecolor = 'black', facecolor = 'orange')

        # annotate ion labels for plot
        log_yrange = np.log10(ylims[1]) - np.log10(ylims[0])
        ax[row][col].annotate(ion, xy=(75, np.power(10, np.log10(ylims[0]) + 0.85*log_yrange)), fontsize=20)
    
        if row == 1:
            ax[row][col].set_xlabel('Impact Parameter (kpc)')
        if col == 0:
            ax[row][col].set_ylabel('Ion Column Density ($\mathrm{cm}^{-2}$)')
        fig.tight_layout()
        plt.savefig('../plots/simulated_column_densities_%s.png'%(sim))


ion_list = ['H I', 'C II', 'C IV', 'Si II', 'N V', 'Si III', 'Si IV', 'Mg II']

sim  = sys.argv[1]

plot_multipanel(ion_list, sim)
