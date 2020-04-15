import matplotlib.pylab as plt
import seaborn as sns
import h5py as h5
import numpy as np



def load_data(property_list, fn = None, ovi_label = None,
                  ion = None, model = None, ray_id = None, redshift = None, flag = None):
    if fn == None:
        fn = '../../data/analyzed_spectra/combined_spectra.h5'

    data = h5.File(fn, 'r')
    
    # always True
    mask = data['impact'][()]  > -99
    if redshift:
        mask = mask & (data[('redshift')][()]      == redshift)
    if model:
        mask = mask & (data[('model')][()]  ==    model)
    if ray_id:
        mask = mask & (data['ray_id'  ][()]  ==   ray_id)
    if ion:
        mask = mask & (data['ion'     ][()]  ==      ion)
    if ovi_label:
        mask = mask & (data['label'   ][()]  == ovi_label)
    if flag:
        mask = mask & (daga['flag'    ][()]  ==      flag)
    return_array = []
    for prop in property_list:
        return_array.append(np.array((data[prop][()] )[mask]))

    data.close()
    return np.array(return_array)


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
        ylims  = (-225, 225)
        ylabel = 'Velocity Offset (km / s)'
    elif yfield == 'bval':
        ylims  = (0.1, 79)
        ylabel = 'Doppler Parameter'
        
    return xlims, ylims, xlabel, ylabel

def plot_data_scatter(ion, xfield = 'impact', yfield = 'col', model = None, ovi_label = None,\
                          redshift = None, ax = None, color = 'black', fn =  None, \
                          label = None, annotate_ion = True, axis_labels = True, marker_size = 20, \
                          marker_style = 's', set_ylim = True):
    if fn == None:
        fn = '../../data/analyzed_spectra/combined_spectra.h5'

    xerrfield = '%s_err'%(xfield)
    yerrfield = '%s_err'%(yfield)
    
    if xfield == 'impact':
        xerrfield = 'impact'

    flagfield = 'flag'
    
    field_list = [xfield, xerrfield, yfield, yerrfield, flagfield]

    xscatter, xerr, yscatter, yerr, flagscatter = load_data(field_list, ion = ion, fn = fn, \
                                      model = model, redshift=redshift, ovi_label = ovi_label)
    
    
    
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
        yscatter = 10**yscatter
        yerr = 10**yerr
        ax.set_yscale('log')
        
        sat = flagscatter == 2
        ax.scatter(xscatter[sat], yscatter[sat],  marker = "^", edgecolor = 'black', facecolor = color, s = marker_size)
        nondetect = flagscatter == 3
        ax.scatter(xscatter[nondetect], yscatter[nondetect],  marker = "v", edgecolor = 'black', facecolor = color, s = marker_size)
        detect = flagscatter == 1
        ax.scatter(xscatter[detect], yscatter[detect],  marker = marker_style, edgecolor = 'black', facecolor = color, label = label, s = marker_size)
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
                                model = None, redshift = None, ovi_label = None, marker_size = 20, color = 'black', \
                                label = None, annotate_ion = True, compare = None, set_ylim = True, fn = None):

    if fn is None:
        fn = '../../data/analyzed_spectra/combined_spectra.h5'
    ncols = int((len(ion_list) + nrows - 1)/ nrows)
    if ax is None:
        fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize=(4*ncols, 3.8*nrows), sharex = True, sharey = False)
    xlims, ylims, xlabel, ylabel = plot_details(xfield, yfield)
#    if nrows == 1: 
#        ax = [ax, ax]

    if compare == 'model':
        models = ['P0', 'tempest']
        ovi_labels = [None, None]
        colors = ['amber', 'windows blue']
        colors = sns.xkcd_palette(colors)
        labels = ['Patient 0', 'Tempest']
        marker_styles = ['s', 'o']

    elif compare == 'ovi':
        if model is None:
            model = 'P0'
        models = [model, model,  model]
        ovi_labels = ['nolow', 'broad', 'narrow']
        colors = ['purple', 'green', 'orange']
        labels = ['No-low', 'Broad', 'Narrow']
        marker_styles = ['s', 's', 's']
    else:
        models = [model]
        colors = [color]
        labels = [label]
        ovi_labels = [ovi_label]
        marker_styles = ['s']
    
    for model, ovi_label, label, color, marker_style in zip (models, ovi_labels, labels, colors, marker_styles):
        for i, ion in enumerate(ion_list):
            ion = ion.replace(" ", "")
            row = int(i / ncols)
            col = i - row*ncols
            
            if nrows == 1: 
                if ncols == 1:
                    figax = ax
                else:
                    figax = ax[col]
            else:
                figax = ax[row][col]
            

            plot_data_scatter(ion, xfield = xfield, yfield = yfield, ax = figax, ovi_label = ovi_label, \
                              model = model, redshift = redshift, fn = fn, \
                              color = color, label = label, marker_size = marker_size, marker_style = marker_style,\
                              annotate_ion = annotate_ion, axis_labels = False, set_ylim = set_ylim)
            
            if row == nrows - 1 and ncols > 1:
                ax[row][col].set_xlabel(xlabel)
            if col == 0 and ncols > 1:
                ax[row][col].set_ylabel(ylabel)

    
    fig.tight_layout()
    return fig, ax    
        
