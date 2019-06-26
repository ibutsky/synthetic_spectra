import matplotlib.pylab as plt
import seaborn as sns
import h5py as h5
import numpy as np



def load_data(property_list, fn = None, use_filtered = True,
                  ion = None, orientation = None, model = None, impact = None, time = None):

    if fn == None:
        if use_filtered:
            fn = '../../data/analyzed_spectra/filtered_spectra.h5'
        else:
            fn = '../../data/analyzed_spectra/combined_spectra.h5'

    data = h5.File(fn, 'r')
    
  #  mask = len(data['impact'].value)* [True]
    mask = data['impact'].value > -99
    if orientation:
        mask = mask & (data['orientation'].value == orientation)
    if model:
        mask = mask & (data['model'      ].value ==       model)
    if impact:
        mask = mask & (data['impact'     ].value ==      impact)
    if time:
        mask = mask & (data['time'       ].value ==        time)
    if ion:
        mask = mask & (data['ion'        ].value ==         ion)

    return_array = []
    for prop in property_list:
        return_array.append(np.array((data[prop].value)[mask]))

    data.close()
    return np.array(return_array)


def apply_column_limits(flag, cols, colerrs, lims, limerrs):
    mask = flag > 8  # check this value
    cols[mask] = lims[mask]
    print(len(cols[mask]))
    colerrs[mask] = np.zeros(len(cols[mask]))
    return cols, colerrs

def plot_details(xfield, yfield):
    if xfield == 'impact':
        xlims  = (0.1, 89)
        xlabel = 'Impact Parameter (kpc)'
    elif xfield == 'bval':
        xlims  = (0.1, 79)
        xlabel = 'Doppler Parameter'

    if yfield == 'col':
        ylims = (1e12, 1e16)
        ylabel = 'Ion Column Density (cm$^{-2}$)'
    elif yfield == 'vel':
        ylims  = (-59, 59)
        ylabel = 'Velocity Offset (km / s)'
    elif yfield == 'bval':
        ylims  = (0.1, 79)
        ylabel = 'Doppler Parameter'
        
    return xlims, ylims, xlabel, ylabel

def plot_data_scatter(ion, xfield = 'impact', yfield = 'col', model = None, orientation = None, \
                           use_filtered = True, time = None, impact = None, ax = None, color = 'black', \
                          label = None, annotate_ion = True, axis_labels = True, marker_size = 20, \
                          set_ylim = True):

    xerrfield = '%serr'%(xfield)
    yerrfield = '%serr'%(yfield)
    
    if xfield == 'impact':
        xerrfield = 'impact'
    if use_filtered:
        flagfield = 'flag'
    else:
        flagfield = 'flag_aodm'
    field_list = [xfield, xerrfield, yfield, yerrfield, flagfield]

    xscatter, xerr, yscatter, yerr, flagscatter = load_data(field_list, ion = ion, use_filtered = use_filtered,\
                                      model = model, orientation = orientation, impact = impact, time = time)
    
    
    # If using filtered data, only keep data points that are measured (above -9999)
    if use_filtered:
        mask = (xscatter > -9999.) & (yscatter > -9999)
        xscatter = xscatter[mask]
        xerr     = xerr[mask]
        yscatter = yscatter[mask]
        yerr     = yerr[mask]
        flagscatter = flagscatter[mask]


    if ax == None:
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(6, 6))
    xlims, ylims, xlabel, ylabel = plot_details(xfield, yfield)
    if axis_labels:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    ax.set_xlim(xlims)
    if set_ylim:
        ax.set_ylim(ylims)

    
    if yfield == 'col':
        if not use_filtered:
            yscatter_lim, yscatter_limerr  = load_data(['col_aodm', 'col_aodm_err'], ion = ion, model = model, \
                                     orientation = orientation, impact = impact, time = time, use_filtered = use_filtered)
            mask = (flagscatter  != 1)
            yscatter[mask] = yscatter_lim[mask]
            yerr[mask]     = yscatter_limerr[mask]

        yscatter = 10**yscatter
        yerr = 10**yerr
        ax.set_yscale('log')
        
        sat = flagscatter == 9
        ax.scatter(xscatter[sat], yscatter[sat],  marker = "^", edgecolor = color, facecolor = 'white', s = marker_size)
        nondetect = flagscatter == 5
        ax.scatter(xscatter[nondetect], yscatter[nondetect],  marker = "v", edgecolor = color, facecolor = 'white', s = marker_size)
        detect = flagscatter == 1
        ax.scatter(xscatter[detect], yscatter[detect],  marker = "s", edgecolor = color, facecolor = 'white', label = label, s = marker_size)
        ax.errorbar(xscatter[detect], yscatter[detect], yerr = yerr[detect], color = 'black', linestyle = '')

        
    else:
        ax.scatter(xscatter, yscatter, edgecolor = color, label = label, facecolor = 'white', s = marker_size)
    ax.errorbar(xscatter, yscatter, yerr = yerr, color = color, facecolor = None, linestyle = '')
    if annotate_ion:
        annotate_ion_name(ax, ion)


def annotate_ion_name(ax, ion_name, fontsize = 18, x_factor = 0.85, y_factor = 0.88, color = 'black'):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ion = ion_name.replace(' ', '')
    if ion  in ['OVI', 'NIII', 'SiIII', 'CIII']:
        x_factor *= 0.95
    elif ion in ['SiIV']:
        x_factor *= 0.92

    if ax.get_xscale() == 'log':
        log_xrange = np.log10(xlim[1]) - np.log10(xlim[0])
        x = np.power(10, np.log10(xlim[0]) + x_factor*log_xrange)
    else:
        x_range = xlim[1] - xlim[0]
        x = xlim[0] + x_factor * x_range

    if ax.get_yscale() == 'log':
        log_yrange = np.log10(ylim[1]) - np.log10(ylim[0])
        y = np.power(10, np.log10(ylim[0]) + y_factor*log_yrange)
    else:
        y_range = ylim[1] - ylim[0]
        y = ylim[0] + y_factor * y_range

    ax.annotate(ion_name, xy = (x,y), fontsize = fontsize, color = color)




    
def plot_multipanel_scatter(ion_list, xfield = 'impact', yfield = 'col', nrows = 2, fig = None, ax = None, \
                                model = None, orientation = None, time = None, impact = None, marker_size = 20,\
                                color = 'black', label = None, annotate_ion = True, compare = None, use_filtered = True, 
                            set_ylim = True):
    ncols = int((len(ion_list) + nrows - 1)/ nrows)
    if ax is None:
        fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize=(4*ncols, 3.8*nrows), sharex = True, sharey = False)
    xlims, ylims, xlabel, ylabel = plot_details(xfield, yfield)
    if nrows == 1: 
        ax = [ax, ax]

    if compare == 'model':
        models = ['anisd', 'stream']
        orientations = [None, None]
        colors = ['amber', 'windows blue']
        colors = sns.xkcd_palette(colors)
        labels = ['Diffusion', 'Streaming']

    elif compare == 'orientation':
        models = [None, None, None]
        orientations = ['face', 'edge_theta0', 'edge_theta1.0']
        colors = ['pale red', 'medium green', 'denim blue']
        colors = sns.xkcd_palette(colors)
        labels = ['face-on', 'theta = 0', 'theta = 1']

    else:
        models = [model]
        orientations = [orientation]
        colors = [color]
        labels = [label]
    
    for orientation, model, label, color in zip (orientations, models, labels, colors):
        for i, ion in enumerate(ion_list):
            ion = ion.replace(" ", "")
            row = int(i / ncols)
            col = i - row*ncols
        
        
            plot_data_scatter(ion, xfield = xfield, yfield = yfield, ax = ax[row][col], \
                              model = model, orientation = orientation, time = time, use_filtered = use_filtered,\
                              impact = impact, color = color, label = label, marker_size = marker_size,\
                              annotate_ion = annotate_ion, axis_labels = False, set_ylim = set_ylim)
            
            if row == nrows - 1:
                ax[row][col].set_xlabel(xlabel)
            if col == 0:
                ax[row][col].set_ylabel(ylabel)

    
    fig.tight_layout()
    return fig, ax    
        
