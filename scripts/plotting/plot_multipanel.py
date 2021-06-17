import matplotlib
#matplotlib.use('Agg')
import sys

import yt


import h5py as h5
import numpy as np
import palettable

from yt.visualization.base_plot_types import get_multi_plot
import matplotlib.colorbar as cb
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LinearSegmentedColormap

#import seaborn as sns

import plotting_tools as pt

#model  = sys.argv[1]
output = 3195

#params = {"text.color" : "white",
#          "xtick.color" : "white",
#          "ytick.color" : "white"}
params = {"xtick.direction": "in",
        "ytick.direction":"in"
}
plt.rcParams.update(params)

field_list = [('gas', 'density'), ('Gas', 'Temperature'), ('Gas', 'relative_velocity_z'), \
              ('gas', 'O_p5_number_density'),\
              ('gas', 'Si_p2_number_density'), ('gas', 'H_p0_number_density')]

#Blues_r
white = (1, 1, 1)
black = (0, 0, 0)
cmap_list = ['magma', 'afmhot', 'GRN','bone',  'dusk', 'purple_mm']
# save for cr eta
# palettable.scientific.diverging.Vik_20.mpl_colormap,
# making custom tempo
tempo = palettable.cmocean.sequential.Tempo_6_r.mpl_colors[1:]
tokyo = palettable.scientific.sequential.Tokyo_20.mpl_colors
tokyo = np.vstack((black, tokyo))
BYW = np.vstack((white, palettable.cartocolors.sequential.BrwnYl_7.mpl_colors))
mint = np.vstack((white, palettable.cartocolors.sequential.Mint_7.mpl_colors))
PBW = np.vstack((white, palettable.colorbrewer.sequential.PuBu_9.mpl_colors))



BYW = palettable.cartocolors.sequential.BrwnYl_7_r.mpl_colors#[1:]
mint = palettable.cartocolors.sequential.Mint_7_r.mpl_colors#[1:]
PBW = palettable.colorbrewer.sequential.PuBu_9_r.mpl_colors#[1:]

cmap_list = ['magma',#LinearSegmentedColormap.from_list('mytempo', tempo),
             palettable.scientific.sequential.LaJolla_20_r.mpl_colormap,
           #  palettable.cartocolors.diverging.Geyser_7.mpl_colormap,
#             palettable.scientific.sequential.Turku_20.mpl_colormap,
                    pt.get_cmap('radial_velocity'),
          #   LinearSegmentedColormap.from_list('Tokyo', tokyo),
          #   palettable.scientific.sequential.Oslo_20.mpl_colormap,
             LinearSegmentedColormap.from_list('ovi', BYW),
             LinearSegmentedColormap.from_list('si', mint),
             LinearSegmentedColormap.from_list('hi', PBW),

             #palettable.scientific.sequential.LaPaz_20_r.mpl_colormap,
            # palettable.scientific.sequential.Davos_20_r.mpl_colormap,
             ]
#zlim_list = [(1e-30, 1e-25), (1e5, 1e8), (5e-3, 5), (1e-21, 1e-15), (1e13, 1e15), (3e12, 1e17)] 
zlim_list = [(1e-28, 5e-25), (3e4, 1.5e6), (-100, 100),  #(1e-2, 1),
            (1e13, 3.2e16), (3.2e12, 1e16), (1e13, 1e20)]

cbar_title_list =[r'$\mathrm{Density}\ (\mathrm{g\ cm^{-3}})$', \
                  r'$\mathrm{Temperature}\ (\mathrm{K})$', \
#                  r'$\mathrm{Metallicity\ } (\mathrm{Z_{\odot}})$',
                    #r'$\mathrm{Radial\ Velocity}\ (\mathrm{km/s})$',
                  r'$\mathrm{LOS\ Velocity}\ (\mathrm{km/s})$',
                  r'$\mathrm{O\ VI\ Column\ Density}\ (\mathrm{cm^{-2}})$', \
                  r'$\mathrm{Si\ III\ Column\ Density}\ (\mathrm{cm^{-2}})$', \
                  r'$\mathrm{H\ I\ Column\ Density}\ (\mathrm{cm^{-2}})$']


# load in simulation data and add ion fields




for i, plot_type in enumerate(['sim', 'col']):
    #initiate figure and axes
    color = 'black'
  #  if i == 1:
  #      plt.rcParams.update(params)
  #      color = 'white'
    orient = 'horizontal'
    model_list = ['P0', 'P0_agncr']
    nrows = len(model_list)
    ncols = 3
    fig, axes, colorbars = get_multi_plot(ncols, nrows, colorbar=None, bw = 4)
    #fig, axes = plt.subplots(nrows = 2, ncols = 3, figsize = (12, 8))
    for row, model in enumerate(model_list):
        frb = h5.File('../../data/simulation_data/multipanel_data_%s_%06d_100_kpc'%(model, output), 'r')
        for col in range(ncols):
            idx = i*3 + col
            field = field_list[idx][1]
            img_data = np.array(frb['%s_z'%field])
            if field.__contains__('velocity'):
                img_data /= 1e5 # converting to km/s
                im = axes[row][col].imshow(img_data, origin = 'lower', vmin = zlim_list[idx][0], vmax = zlim_list[idx][1])
            else:
                    im = axes[row][col].imshow(img_data, origin = 'lower', norm = LogNorm(),\
                               vmin = zlim_list[idx][0], vmax = zlim_list[idx][1])
            im.set_cmap(pt.get_cmap(field, listed = False))#cmap_list[idx])
            if field == ('density'):
                cmap = palettable.cmocean.sequential.Tempo_6_r.mpl_colormap
                #cmap = palettable.scientific.sequential.Turku_8.mpl_colormap
                #cmap = palettable.scientific.sequential.LaPaz_8.mpl_colormap
                #cmap = palettable.scientific.sequential.Davos_8.mpl_colormap
                #cmap = palettable.cmocean.sequential.Deep_6_r.mpl_colormap
                #cmap = palettable.cmocean.sequential.Matter_20_r.mpl_colormap


                im.set_cmap(cmap)
           # elif field.__contains__('Temperature'):
            #    im.set_cmap('afmhot')

            axes[row][col].set_xticks([0, 200, 400, 600, 800])
            axes[row][col].set_yticks([0, 200, 400, 600, 800])

            axes[row][col].set_xticklabels('')
            axes[row][col].set_yticklabels('')
            axes[row][col].xaxis.set_visible(True)
            axes[row][col].yaxis.set_visible(True)
            axes[row]

            cbax = inset_axes(axes[row][col], width = "90%", height = "3%", loc = 9)
            cbar = fig.colorbar(im, cax=cbax, orientation = 'horizontal')
            cbar.set_label(cbar_title_list[idx], color = color)
   # axes[0][0]
    axes[0][0].text(50, 50, 'P0', color = 'black', bbox=dict(facecolor = 'white', edgecolor = 'black', boxstyle='round, pad=0.5' ))
    axes[1][0].text(50, 50, 'P0 + CR', color = 'black', bbox=dict(facecolor = 'white', edgecolor = 'black', boxstyle='round, pad=0.5' ))
    
    axes[0][0].plot([550, 750], [100, 100], linewidth = 4, color = 'black')
    axes[0][0].text(600, 130, '25 kpc', color = 'black')
    # note: resolution is 800 x 800 which spans -50 to 50 kpc in both directions

    # And now we're done!
    #fig.tight_layout(pad=2)
   # fig.tight_layout(rect=[0,0,.8,1])

    fig.savefig("../../plots/multipanel_%06d_%s_z.png"%(output, plot_type), dpi = 300)
