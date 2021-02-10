import matplotlib
#matplotlib.use('Agg')
import sys

import yt


import h5py as h5
import numpy as np
import palettable
import seaborn as sns


from yt.visualization.base_plot_types import get_multi_plot
import matplotlib.colorbar as cb
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LinearSegmentedColormap, ListedColormap


#model  = sys.argv[1]
output = 3195

params = {"text.color" : "white",
          "xtick.color" : "white",
          "ytick.color" : "white"}
#plt.rcParams.update(params)

field_list = [('gas', 'density'), ('Gas', 'Temperature'), ('Gas', 'metallicity2'), \
              ('gas', 'O_p5_number_density'),\
              ('gas', 'Si_p2_number_density'), ('gas', 'H_p0_number_density')]

#Blues_r
white = (1, 1, 1)
black = (0, 0, 0)
cmap_list = ['magma', 'afmhot', 'GRN','bone',  'dusk', 'purple_mm']
# save for cr eta
# palettable.scientific.diverging.Vik_20.mpl_colormap,
# making custom tempo
tempo = palettable.cmocean.sequential.Tempo_6.mpl_colors[:]
tokyo = palettable.scientific.sequential.Tokyo_20.mpl_colors
tokyo = np.vstack((black, tokyo))
BYW = np.vstack((white, palettable.cartocolors.sequential.BrwnYl_7.mpl_colors))
mint = np.vstack((white, palettable.cartocolors.sequential.Mint_7.mpl_colors))
PBW = np.vstack((white, palettable.colorbrewer.sequential.PuBu_9.mpl_colors))

cmap_ovi = ListedColormap(palettable.cartocolors.sequential.BrwnYl_7.mpl_colors)
cmap_siii = ListedColormap(palettable.cartocolors.sequential.Mint_7.mpl_colors)
cmap_hi  = ListedColormap(sns.color_palette("Blues"))

BYW = palettable.cartocolors.sequential.BrwnYl_7_r.mpl_colors#[1:]
mint = palettable.cartocolors.sequential.Mint_7_r.mpl_colors#[1:]
PBW = palettable.colorbrewer.sequential.PuBu_9_r.mpl_colors#[1:]

cmap_list = [palettable.cmocean.sequential.Tempo_6.mpl_colormap,
             palettable.scientific.sequential.LaJolla_20_r.mpl_colormap,
             palettable.cartocolors.diverging.Geyser_7.mpl_colormap,
             cmap_ovi, cmap_siii, cmap_hi
             #palettable.scientific.sequential.Turku_20.mpl_colormap,
             #LinearSegmentedColormap.from_list('Tokyo', tokyo),
             #palettable.scientific.sequential.Oslo_20.mpl_colormap,
             
            # LinearSegmentedColormap.from_list('ovi', BYW),
            # LinearSegmentedColormap.from_list('si', mint),
             #LinearSegmentedColormap.from_list('hi', PBW),

             #palettable.scientific.sequential.LaPaz_20_r.mpl_colormap,
            # palettable.scientific.sequential.Davos_20_r.mpl_colormap,
             ]
#zlim_list = [(1e-30, 1e-25), (1e5, 1e8), (5e-3, 5), (1e-21, 1e-15), (1e13, 1e15), (3e12, 1e17)] 
zlim_list = [(3e-29, 1e-25), (3e4, 1e6), (1e-2, 1), (1e13, 3e15), (1e13, 3e16), (1e13, 1e21)]

cbar_title_list =[r'$\mathrm{Density}\ (\mathrm{g\ cm^{-3}})$', \
                  r'$\mathrm{Temperature}\ (\mathrm{K})$', \
                  r'$\mathrm{Metallicity\ } (\mathrm{Z_{\odot}})$',\
                  r'$\mathrm{O\ VI\ Column\ Density}\ (\mathrm{cm^{-2}})$', \
                  r'$\mathrm{Si\ III\ Column\ Density}\ (\mathrm{cm^{-2}})$', \
                  r'$\mathrm{H\ I\ Column\ Density}\ (\mathrm{cm^{-2}})$']


# load in simulation data and add ion fields




for i, plot_type in enumerate(['sim', 'col']):
    #initiate figure and axes
    color = 'black'
   # if i == 1:
   ##     plt.rcParams.update(params)
     #   color = 'white'
    orient = 'horizontal'
    model_list = ['P0', 'P0_agncr']
    nrows = len(model_list)
    ncols = 3
    fig, axes, colorbars = get_multi_plot(ncols, nrows, colorbar=None, bw = 4)
    for row, model in enumerate(model_list):
        frb = h5.File('../../data/simulation_data/multipanel_data_%s_%06d'%(model, output), 'r')
        for col in range(ncols):
            idx = i*3 + col
            img_data = np.array(frb[field_list[idx][1]])
            im = axes[row][col].imshow(img_data, origin = 'lower', norm = LogNorm(),\
                               vmin = zlim_list[idx][0], vmax = zlim_list[idx][1])
            im.set_cmap(cmap_list[idx])

            axes[row][col].xaxis.set_visible(False)
            axes[row][col].yaxis.set_visible(False)

            cbax = inset_axes(axes[row][col], width = "90%", height = "3%", loc = 9)
            cbar = fig.colorbar(im, cax=cbax, orientation = 'horizontal')
            cbar.set_label(cbar_title_list[idx], color = color)
            
    axes[0][0].text(50, 50, 'P0', color = 'black', bbox=dict(facecolor = 'white', edgecolor = 'black', boxstyle='round, pad=0.5' ))
    axes[1][0].text(50, 50, 'P0 + CR', color = 'black', bbox=dict(facecolor = 'white', edgecolor = 'black', boxstyle='round, pad=0.5' ))

    # And now we're done!
    fig.savefig("../../plots/multipanel_%06d_%s.png"%(output, plot_type), dpi = 300)
