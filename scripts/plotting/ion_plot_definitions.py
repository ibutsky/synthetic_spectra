import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import numpy as np
import h5py as h5
import sys
import os.path

def return_ylims(ion):
    if ion == 'H I':
        ylims = (1e11, 1e21)
    elif ion == 'O VI':
        ylims = (1e10, 1e17)
    elif ion == 'Mg II':
        ylims = (1e2, 1e17)
    elif ion == 'C II':
        ylims = (1e5, 1e19)
    elif ion == 'C III':
        ylims = (1e7, 1e17)
    elif ion == 'C IV':
        ylims = (1e7, 1e16)
    elif ion == 'Si II':
        ylims = (1e2, 1e17)
    elif ion == 'Si III':
        ylims = (1e5, 1e16)
    elif ion == 'Si IV':
        ylims = (1e5, 1e16)
    else:
        ylims = (1e5, 1e20)
    return ylims 

def return_observational_threshold(ion):
    ion = ion.replace(" ", "")
    if ion == 'HI':
        return np.power(10, 13)
    elif ion == 'O VI': 
        return np.power(10, 13.5)
    elif ion == 'Si II': 
        return np.power(10, 13) 
    elif ion == 'Si III': 
        return np.power(10, 12.8)
    elif ion == 'Si IV': 
        return np.power(10, 13.3) 
    elif ion == 'C II': 
        return np.power(10, 13.3)
    elif ion == 'C III':
        return np.power(10, 13)
    else:
        return np.power(10, 13)


def return_ion_prefix(ion):
    if ion == 'H I':
        field = 'H_p0'
    elif ion == 'O VI':
        field = 'O_p5'
    elif ion == 'Si II':
        field = 'Si_p1'
    elif ion == 'Si III':
        field = 'Si_p2'
    elif ion == 'Si IV':
        field = 'Si_p3'
    elif ion == 'C II':
        field = 'C_p1'
    elif ion == 'C III':
        field = 'C_p2'
    elif ion == 'C IV':
        field = 'C_p3'
    elif ion == 'Mg II':
        field = 'Mg_p1'
    elif ion == 'N V':
        field = 'N_p4'
    else:
        print(ion)
    return field

def return_field_name(ion, field, full_name = True):
    ion_prefix = return_ion_prefix(ion)
    field_name = "%s_%s"%(ion_prefix, field)
    if full_name:
        field_name = ('gas', field_name)
    return field_name

def generate_ion_field_list(ion_list, field, full_name = True):
    field_name_list = []
    for ion in ion_list:
        field_name_list.append(return_field_name(ion, field, full_name))
    return field_name_list
                        

def make_profiles(r_arr, cdens_arr, r_bins, n_bins):
    bin_ids = np.digitize(r_arr, r_bins)
    median = np.zeros(len(r_bins))
    std = np.zeros(len(r_bins))
    for i in np.arange(n_bins):
        bin_id = i + 1
        sample = cdens_arr[bin_ids == bin_id]
        median[i] = np.median(sample)
        std[i] = np.std(sample)
        
    return median, std

def make_covering_fraction_profiles(r_arr, cdens_arr, r_bins, n_bins, threshold):
    bin_ids = np.digitize(r_arr, r_bins)
    covering_fraction_profile_data = np.zeros(len(r_bins))
    for i in np.arange(n_bins):
        bin_id = i + 1
        sample = cdens_arr[bin_ids == bin_id]
        above_threshold = sample[sample > threshold]
        covering_fraction_profile_data[i] = len(above_threshold) / len(sample)
    return covering_fraction_profile_data

def make_median_and_cfrac_profiles(r_arr, cdens_arr, r_bins, n_bins, threshold):
    bin_ids = np.digitize(r_arr, r_bins)
    median = np.zeros(len(r_bins))
    std = np.zeros(len(r_bins))
    covering_fraction_profile_data = np.zeros(len(r_bins))
    for i in np.arange(n_bins):
        bin_id = i + 1
        sample = cdens_arr[bin_ids == bin_id]
        median[i] = np.median(sample)
        std[i] = np.std(sample)
        above_threshold = sample[sample > threshold]
        if (len(sample)) > 0:
            covering_fraction_profile_data[i] = len(above_threshold) / len(sample)
        else:
            covering_fraction_profile_data[i] = 0
    return median, std, covering_fraction_profile_data


def plot_hist2d(ax, r_arr, cdens_arr, rmax,  ylims, cmap='GnBu', vmin=1, vmax=None):
    nbins = 400
    xbins = np.linspace(0, rmax, nbins)
    ybins = 10**np.linspace(np.log10(ylims[0]), np.log10(ylims[1]), nbins)
    counts, x_edge, y_edge = np.histogram2d(r_arr, cdens_arr, bins=(xbins, ybins))
    x_bin_center = ((x_edge[1:] + x_edge[:-1]) / 2).reshape(nbins-1,1)
    # normalize counts in x-space to remove out linear increase in counts with 
    # radius due to circles of constant impact parameter
    counts /= x_bin_center 
    ax.set_yscale('log')
    #im = ax.pcolormesh(xbins, ybins, counts.T, vmin=counts.min(), vmax=counts.max(), cmap='magma', norm=LogNorm())
    print(counts.max())
    if vmax == None:
        vmax = counts.max()
    im = ax.pcolormesh(xbins, ybins, counts.T, vmin=vmin, vmax=vmax, cmap=cmap, norm=LogNorm())
    return im



def load_r_cdens(fname, ion, underscore = False):
    r_arr = []
    cdens_arr = []

    frb = h5.File(fname, 'r')
    for axis in ['x', 'y', 'z']:
        cname = "%s_%s"%(ion.replace(" ", ""), axis)
           
        if cname in frb.keys():
            r_arr = np.concatenate((r_arr, frb['radius'].value))
            cdens_arr = np.concatenate((cdens_arr, frb[cname].value))
    return r_arr, cdens_arr


def make_projection(ds, axis, ion_fields, center, width):
    p = ds.proj(ion_fields, axis, weight_field=None, center=center, method='integrate')
    return p.to_frb(width, 800, center=center)


def generate_mask(data, ion_list, model =  None, orientation = None):
    mask = np.zeros(len(data['impact']), dtype=bool)
    if ion_list != 'all':
        for ion in ion_list:
            ion = ion.replace(" ", "")
            mask = mask | ((data['ion'].value == ion))
    if model: 
        mask = mask & (data['model'].value == model)
    if orientation:
        mask = mask & (data['orientation'].value == orientation)
        
    return mask


def generate_masked_data(ion_list, xfield, yfield, flagfield, model = None, orientation = None):
    mask = generate_mask(ion_list, model = model, orientation = orientation)
    if xfield == 'impact':
        xlist = impact_list
        xerr  = np.zeros(len(impact_list)) 
    elif xfield == 'bval':
        xlist = bval_list
        xerr = sigbval_list
    elif xfield == 'vel':
        xlist = vel_list
        xerr = sigvel_list
    if xfield == 'col':
        xlist = col_list
        xerr = sigcol_list
        
    if yfield == 'col':
        ylist = col_list
        yerr = sigcol_list
    elif yfield == 'vel':
        ylist = vel_list
        yerr = sigvel_list
    elif yfield == 'impact':
        ylist = impact_list
        yerr  = np.zeros(len(impact_list)) 
    elif yfield == 'bval':
        ylist = bval_list
        yerr = sigbval_list
        
    return xlist[mask], xerr[mask], ylist[mask], yerr[mask], flagfield[mask]

def plot_details(xfield, yfield):
    if xfield == 'impact':
        xlims  = (0.1, 69)
        xlabel = 'Impact Parameter (kpc)'
    elif xfield == 'bval':
        xlims  = (0.1, 79)
        xlabel = 'Doppler Parameter'

    if yfield == 'col':
        #ylims  = (10, 17)
        ylims = (1e12, 1e15)
        ylabel = 'Column Density (cm$^{-2}$)'
    elif yfield == 'vel':
        ylims  = (-59, 59)
        ylabel = 'Velocity Offset (km / s)'
        
    return xlims, ylims, xlabel, ylabel

def apply_column_limits(flag, cols, colerrs, lims, limerrs):
    mask = flag > 8  # check this value
    cols[mask] = lims[mask]
    print(len(cols[mask]))
    colerrs[mask] = np.zeros(len(cols[mask]))
    return cols, colerrs


def load_veeper_ion_column(ion, model, orientation = None):
    data = h5.File('../data/combined_spectra.h5', 'r')
    flag_list = data['flag'].value
    impact_list = data['impact'].value
    col_list = data['col'].value
    sigcol_list = data['colerr'].value
    col_list, sigcol_list = apply_column_limits(flag_list, col_list,\
                                            sigcol_list, data['col_lim'], data['col_limerr'])

    mask = generate_mask(data, [ion], model = model, orientation = orientation)
    return impact_list[mask], col_list[mask], sigcol_list[mask], flag_list[mask]

