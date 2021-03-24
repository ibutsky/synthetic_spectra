import matplotlib.pylab as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

from yt import YTArray
import seaborn as sns
import palettable
import h5py as h5
import numpy as np
import csv
from astropy import constants as const
c = const.c.to('km/s').value



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

def load_csv_data(csvfile = None):
    if csvfile is None:
        csvfile = '/Users/irynabutsky/Downloads/merged_solutions_with_upper_limits.csv'
    data = csv.reader(open(csvfile, 'rt'))

    ion_list, logN_list, logNerr_list, b_list, berr_list, vcen_list, \
        vcen_sd_list, ray_id_list, model_list, detection_list = [], [],[],[],[],[],[],[],[],[]

    counter = 0
    dummy = -9999
    for row in data:
        if counter > 0:
            detection = row[9]
            detection_list.append(detection)
            ion_list.append(row[1])
            logN_list.append(float(row[2]))
            model_list.append(row[6][:-4])
            ray_id_list.append(int(row[6][-4:]))
            if detection == 'detection':
                logNerr_list.append(float(row[3]))
                b_list.append(float(row[4]))
                berr_list.append(float(row[5]))
                vcen_list.append(float(row[7]))
                vcen_sd_list.append(float(row[8]))
            else:
                logNerr_list.append(dummy)
                b_list.append(dummy)
                berr_list.append(dummy)
                vcen_list.append(dummy)
                vcen_sd_list.append(dummy)
        counter += 1
                          
    ion_list = np.array(ion_list)
    logN_list = np.array(logN_list)
    logNerr_list = np.array(logNerr_list)
    b_list = np.array(b_list)
    berr_list = np.array(berr_list)
    vcen_list = np.array(vcen_list)
    vcen_sd_list = np.array(vcen_sd_list)
    ray_id_list = np.array(ray_id_list)
    model_list = np.array(model_list)
    detection_list = np.array(detection_list)

#    zsys = 0.25
#    vel_list = c * ((1 + zcen_list) / (1 + zsys) - 1)

    # make O VI labels
    ovi_label_list = np.array(len(vcen_list)*[None])
    for i in range(max(ray_id_list)):
        for model in ['P0', 'P0agncr']:
            mask = (ray_id_list == i) & (model_list == model) & (detection_list == 'detection')
            ovi_mask = mask & (ion_list == 'OVI') 
            siIII_mask = mask & (ion_list == 'SiIII') 
            #print(i, model, ion_list[ovi_mask], ion_list[siIII_mask])
            if len(ion_list[ovi_mask]) == 0:
                continue
            if len(ion_list[siIII_mask]) == 0:
                ovi_label_list[ovi_mask] = 'nolow'
            else:
                o_vel_list = vcen_list[ovi_mask]
                si_vel_list = vcen_list[siIII_mask]
                for j, o_vel in enumerate(o_vel_list):
                    dvel_list = o_vel - si_vel_list
                    ovi_label_mask = ovi_mask & (vcen_list == o_vel)
                    if min(np.abs(dvel_list)) > 35:
                        ovi_label_list[ovi_label_mask] = 'nolow'
                    elif b_list[ovi_label_mask] > 30:
                        ovi_label_list[ovi_label_mask] = 'broad'
                    else:
                        ovi_label_list[ovi_label_mask] = 'narrow'
            
    idx, impacts = np.loadtxt('../../data/random_sightlines.dat', unpack = True, skiprows = 1, usecols = (0, 1))
    impact_list = np.array([])
    for i, ray_id in enumerate(ray_id_list):
        mask = idx == ray_id
        r = impacts[mask][0]
        impact_list = np.append(impact_list, r)
        
    return ion_list, logN_list, logNerr_list, b_list, berr_list, vcen_list, \
            vcen_sd_list, ray_id_list, model_list, impact_list, ovi_label_list, detection_list
            
def get_total_column(ray_id_list, log_col_list):
    total_col = np.array(len(ray_id_list)*[0.0])
    for ray_id in range(0, max(ray_id_list)):
        mask = ray_id_list == ray_id
        cols = np.sum(np.power(10, log_col_list[mask]))
        total_col[mask] = np.log10(cols)
    return total_col
    
    
def convert_to_vel(wl, w0):
    return (wl-w0) / w0 * 2.9979e5

def los_vel(vx, vy, vz, bv = [0, 0, 0], normal = [0, 0, 1]):
    vx = np.array(vx - bv[0])
    vy = np.array(vy - bv[1])
    vz = np.array(vz - bv[2])

    normal_mag = np.linalg.norm(normal)
    normal = np.divide(normal, normal_mag)

    v_dot_norm = vx*normal[0] + vy*normal[1] + vz*normal[2]
    return v_dot_norm
    

def load_sightline_scatter_data(sim, ray_id, output = 3195):
    fn = '../../data/unanalyzed_spectra/ray_%s_%i_%i.h5'%(sim, output, ray_id)
    plot_data = h5.File(fn, 'r')['grid']

    l = YTArray(plot_data['l'], 'cm')
    l = np.array(l.in_units('kpc'))

    temperature = np.array(plot_data['temperature'])
    density     = np.array(plot_data['density'])
    metallicity = np.array(plot_data['metallicity'])*77.22007722007721  # converting from code units to zsun
    vx          = np.array(plot_data['relative_velocity_x']) / 1e5 # converting to km/s
    vy          = np.array(plot_data['relative_velocity_y']) / 1e5 # converting to km/s
    vz          = np.array(plot_data['relative_velocity_z']) / 1e5 # converting to km/s

    vlos        = np.array(plot_data['velocity_los']) / 1e5
    dl = np.array(plot_data['dl'])
    
    # O VI and H I column densities
    oden = np.array(plot_data['O_p5_number_density'])
    ocol  = dl * np.array(oden)
    sicol = dl* np.array(plot_data['Si_p2_number_density'])
    hcol  = dl * np.array(plot_data['H_p0_number_density'])
    
    return l, temperature, density, metallicity, vlos, ocol, sicol

def plot_details(xfield, yfield, log=True):
    if xfield == 'impact':
        xlims  = (0.1, 89)
        xlabel = 'Impact Parameter (kpc)'

    elif xfield == 'bval':
        xlims  = (0.1, 79)
        xlabel = 'Doppler Parameter'
        if log:
            xlims = np.log10(xlims)
    if yfield == 'col' or yfield == 'total_col':
        ylims = (1e12, 1e16)
        ylabel = 'Ion Column Density (cm$^{-2}$)'
        if log:
            ylims = (12, 16)
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
    
    #print(xscatter[1])
    xscatter = np.log10(xscatter)
    mask = (xscatter > -9999.) & (yscatter > -9999) & (yerr < 0.5)   #TEMPORARY yerr flag
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

    
    if yfield == 'col' or yfield == 'total_col':
#        yscatter = 10**yscatter
#        yerr = 10**yerr
#        ax.set_yscale('log')
        
        sat = flagscatter == 2
        ax.scatter(xscatter[sat], yscatter[sat],  marker = "^", edgecolor = 'black', facecolor = color, s = marker_size)
        nondetect = flagscatter == 3
        ax.scatter(xscatter[nondetect], yscatter[nondetect],  marker = "v", edgecolor = 'black', facecolor = color, s = marker_size)
        detect = flagscatter == 1
        ax.errorbar(xscatter[detect], yscatter[detect], yerr = yerr[detect], color = 'black', linestyle = '', zorder = 0, alpha = 0.2)
        ax.scatter(xscatter[detect], yscatter[detect],  marker = marker_style, edgecolor = 'black', facecolor = color, label = label, s = marker_size)


        
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
        models = ['P0', 'P0_agncr']
        ovi_labels = [None, None]
        colors = ['amber', 'windows blue']
        colors = sns.xkcd_palette(colors)
        labels = ['Patient 0', 'P0+agncr']
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

def get_cmap(field):
    if field =='density':
        cmap = palettable.cmocean.sequential.Tempo_6.mpl_colormap
    elif field == 'pressure':
        cmap = 'magma'
    elif field == 'temperature' or field == 'Temperature':
        cmap = palettable.scientific.sequential.LaJolla_20_r.mpl_colormap
    elif field == 'cr_eta':
        cmap = palettable.scientific.sequential.Tokyo_20.mpl_colormap
    elif field == 'cr_pressure':
        cmap = palettable.scientific.sequential.Turku_20.mpl_colormap
    elif field.__contains__('velocity') or field == 'radial_velocity':
        cmap = palettable.scientific.diverging.Vik_20.mpl_colormap
    elif field == 'magnetic_field_strength':
        cmap = palettable.scientific.sequential.LaPaz_20.mpl_colormap
    elif field.__contains__('H_p'):
        cmap = ListedColormap(sns.color_palette("Blues", 7))
    elif field.__contains__('Si_p'):
        cmap = ListedColormap(palettable.cartocolors.sequential.Mint_7.mpl_colors)
    elif field.__contains__('O_p'):
        cmap = ListedColormap(palettable.cartocolors.sequential.BrwnYl_7.mpl_colors)
    elif field == 'metallicity' or field == 'metallicity2':
        cmap = palettable.cartocolors.diverging.Geyser_7.mpl_colormap
    else:
        cmap = 'viridis'
        print(field, "doesn't have a designated colormap")
    return cmap

def get_palette(field):
    if field =='density':
        cmap = palettable.cmocean.sequential.Tempo_6.mpl_colors
    elif field == 'temperature' or field == 'Temperature':
        cmap = palettable.scientific.sequential.LaJolla_20_r.mpl_colors
    elif field == 'cr_eta':
        cmap = palettable.scientific.sequential.Tokyo_20.mpl_colors
    elif field == 'cr_pressure':
        cmap = palettable.scientific.sequential.Turku_20.mpl_colors
    elif field == 'velocity_magnitude' or field == 'velocity_los':
        cmap = palettable.scientific.diverging.Vik_20.mpl_colors
    elif field == 'magnetic_field_strength':
        cmap = palettable.scientific.sequential.LaPaz_20.mpl_colors
    elif field.__contains__('H_p') or field == 'H I':
        cmap = sns.color_palette("Blues", 7)
    elif field.__contains__('Si_p') or field == 'Si III':
        cmap = palettable.cartocolors.sequential.Mint_7.mpl_colors
    elif field.__contains__('O_p') or field == 'O VI':
        cmap = palettable.cartocolors.sequential.BrwnYl_7.mpl_colors
    elif field == 'metallicity' or field == 'metallicity2':
        cmap = palettable.cartocolors.diverging.Geyser_7.mpl_colors
    else:
        cmap = palettable.cmocean.sequential.Tempo_6.mpl_colors
        print(field, "doesn't have a designated color palette")
    return cmap
    
#def get_palette(ion):
#    if ion == 'O VI':
 #       palette = palettable.cartocolors.sequential.BrwnYl_7.mpl_colors
  #  elif ion == 'Si III':
  #      palette = palettable.cartocolors.sequential.Mint_7.mpl_colors
  #  elif ion == 'H I':
   #     palette = sns.color_palette("Blues")
  #  return palette

