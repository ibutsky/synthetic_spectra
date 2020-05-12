import matplotlib.pylab as plt
import numpy as np
import scipy
from scipy import interpolate
import sys
import os

from linetools.isgm.abscomponent import AbsComponent
import json
from pyigm.guis import igmguesses
from linetools.spectra.io import readspec


import spectrum_analysis_tools as spa

# inwave = REST wavelength for rest-frame EW; observed wavelength for observed EW
# inspec = flux
# inerror = error
# vrange = [vmin, vmax] velocity range for integration
# w0 = rest wavelength of transition
# f0 = fval of transition
# limit = how many sigma limit do we want, for limits?
# norm = the flux value to normalise the spectrum by
# nstapix -- the number of pixels that need to be below the saturation limit or error before we trigger the saturation flag. 
#
# Keywords
# ctnorm -- normalise the spectra by the median flux in +/-600km/s (excluding the line region)
#
# OUTPUTS
# EW -- 2-element array of eqw and eqwerror. 
# colm -- 2-element array of linear column density and error
#
# Note -- Return 2-element arrays with error included, since 2 keywords eqw and eqwerr are
#         ambiguous to IDL and so not allowed.
#inspec
# REVISION HISTORY
#
def deriv(x):
        # return scipy.misc.derivative(x)
        return np.diff(np.append(x, x[-1]))

def ct_ew2n(ew, wrest, fval):
        # assume optically thin case, and compute column density from EW. See, e.g. Savage & Sembach (1996) eqn 3
        # EW in mA
        # wrest in A
        # fval is from morton

        return 1.13e17*ew/(fval * wrest**2)

def eqwrange(ion, wave, spec, error, vrange, w0, f0, \
              limit = 1., norm = 1.0, ctnorm = 1, nsatpix = 3, \
              silent = 0, plots = 0, sig_limit = 3, sat_limit = 0.1, \
             plot_dir = '.'):


        ### DEFAULT SETTING FOR FLAG. 1 = ok.
        flag = 1

        if norm:
                spec = spec / norm
                error = error / norm
        if ctnorm:
                #grab the median flux value in the +/-600km/s range, outside the line area
                vel = (wave-w0)/w0*2.9979e5
                qq = ((vel > -1200) & (vel < vrange[0])) | ((vel > vrange[1]) & (vel < 1200))
                med = np.median(spec[qq])
                spec  /= med
                error /= med
        else:
                print('warning - assumes data and errors are already continuum normalized')

        # convert wavelengths to velocty
        vv = (wave-w0) / w0 * 2.9979e5
        # find points in range
        # get rid of infinity values down the line
        # note: results sensitive to the value we set spec to
        # sort of converges but not really
        spec[spec <= 0] = 0.01
        iv = (vv > vrange[0]) & (vv < vrange[1]) 
        if len(vv[iv]) == 0:
                # Something is wrong. We were given a spectrum with no points in the velocity range. Pull out!!
                print("Spectrum has no points in the velocity range")
                ew = [-999.9, -999.9]
                colm = [-999.9, -999.9]
                return

        sat = (spec[iv] < sat_limit) | (spec[iv] < error[iv])
        if len(spec[iv][sat]) > 0:
                spec[iv][sat] = error[iv][sat]
                if len(spec[iv][sat]) > nsatpix:
                        flag += 8  #If there are more pixels than we set, then trip the saturated flag

        # now that we've calculated the saturation flag, we can filter spec[iv] 
        # down even more to get more accurate column density estimates
        # iv =  iv & (spec > sat_limit)	

        # plotting the flux vs velocities for all pixels and those used in fit
        if plots:
                fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize=(10, 5), sharex = False, sharey = False)
                ax[0].errorbar(vv, spec, yerr = error)
                ax[0].errorbar(vv[iv], spec[iv], yerr = error[iv])
                ax[0].axhline(y = 0, color = 'red', linestyle = 'dashed')
                ax[0].set_xlim(vrange[0], vrange[1])
                ax[0].set_ylim(-0.1, 1.5)
                ax[0].set_xlabel('Relative Velocity (km/s)')
                ax[0].set_ylabel('Normalized Flux')

        dwave = deriv(wave[iv])
        d_eqw = dwave * (1. - spec[iv])  # * 1000.  # why 1000?
        
        eqw = np.sum(d_eqw)

        # dtermine centroid of equivalent width
        f_eqw = np.cumsum(d_eqw) / np.max(np.cumsum(d_eqw)) #cumulative eqw fraction
        # interpolate for velocity at which fraction of eqw = 0.5
        # we'll take this as the velocity centroid
        # trying something new: 
        vcent = np.interp(0.5, f_eqw, vv[iv])
        eqwerr = np.sqrt(np.sum((deriv(wave[iv])**2  * error[iv])**2) )
        EW = [eqw, eqwerr]

        # Upper limit flag
        if eqw > sig_limit*eqwerr:
                flag += 4
                
        #now convert data points to tau0 -> column density using input line data
        tau = -1. * np.log(spec)
        f_tau = np.cumsum(tau[iv]) / max(np.cumsum(tau[iv])) # cumulative opticaldepth fraction
        tau_cent = np.interp(0.5, f_tau, vv[iv])        
        
        # calculate velocity width as moment of tau dist
        top = np.sum((vv[iv] - tau_cent)**2 * tau[iv] * deriv(vv[iv])) # integral of (v-vbar)^2 * tau_a(v) * dv
        bottom = np.sum(tau[iv] * deriv(vv[iv])) # integral of tau_a(v) * dv
        sigma_v = (top / bottom)**0.5 # this is sigma as defined on pg 692 of heckman et al. 2002
        velcent = tau_cent
        velwidth = np.sqrt(2.) * sigma_v
        
        nv = tau / 2.654e-15 / f0 / w0 # in units cm^-2 / (km s^-1), SS91
                                       # column density per unit velocity
        n = nv * deriv(vv) # column density per bin obtained by multiplying differential Nv by bin width
        # when spec[iv] == 0, tau and therefore n and col = inf
        # fixed this by specifying the [iv] mask to have spec > 0
        for ii in range(len(n[iv])):
                if np.isinf(n[iv][ii]):
                        print(nv[iv][ii], spec[iv][ii])

        tauerr = error / spec
        nerr = tauerr / 2.654e-15 / f0 / w0 * deriv(vv)
        if plots:
                ax[1].errorbar(vv, n,  yerr = nerr)#, drawstyle = 'steps-mid-')
                ax[1].errorbar(vv[iv], n[iv], yerr = nerr[iv])#, drawstyle = 'steps-mid-')
                ax[1].set_yscale('log')
                ax[1].set_ylim(1e10, 1e15)
                ax[1].set_xlim(-250, 250)
                ax[1].set_xlabel('Relative Velocity (km/s)')
                ax[1].set_ylabel('Column Density (cm$^{-2}$)')
                fig.tight_layout()
                plt.savefig('%s/%s_%.1f_flag_%i.png'%(plot_dir, ion, w0, flag))

        col = np.sum(n[iv])
        colerr = np.sqrt(np.sum(nerr[iv]**2))
        colm = [col, colerr]
        
        if not silent:
                if not limit:
                        print('Direct N = %e +/- %e'%(col, colerr))
                        print('Direct N = %e, +/- %e'%(np.log10(col), np.alog10(col+colerr) - alog10(col)))
                        print('W_lambda = %e +/- %e mA over %s km/s'%(eqw, eqwerr, str(vrange)))
                else:
                        print('Limit N = (%e) +/- %e'%(col, colerr*limit))
                        print('Limit N = (%e) +/- %e'%(np.log10(np.abs(col)), np.log10(colerr * limit)))
                        print('W_lambda = (%e) +/- %e mA over %s km/s'%(eqw, eqwerr*limit, str(vrange)))
        return eqw, eqwerr, np.log10(np.abs(col)), np.log10(col+colerr)-np.log10(col), flag, velcent, velwidth


def find_ion_limits(ion, fn, restwave = 0, redshift = 0, sig_limit = 3, \
                    vrange = (-200, 200), silent = 0, plots = 0, plot_dir = '.', sat_limit = 0.1):
    ion_names = np.loadtxt('../../dla.lst', unpack=True, skiprows = 1, usecols=(1), dtype = 'str')
    ion_wls, ion_fs = np.loadtxt('../../dla.lst', unpack=True, skiprows = 1,usecols=(0, 3))
    ion_wls_int = ion_wls.astype(int)
    
    wl, flux, ferr = spa.load_spec_from_fits(fn)
    
    eqw_list = []
    eqwerr_list = []
    col_list = []
    colerr_list = []
    flag_list = []
    velcent_list = []
    velwidth_list = []
    
    ion = ion.replace(" ", "")
    if restwave == 0:
            w0 = spa.restwave(ion, z = 0) # REST wavelength for rest-frame EW; observed wavelength for observed EW
            w_obs = spa.restwave(ion, z = redshift)
    else:
            w0 = restwave
            w_obs = restwave * (redshift + 1)
    ion_mask  = (ion_names == ion)
    wl_mask =  (ion_wls_int[ion_mask] == int(w0))
    f0 = ion_fs[ion_mask][wl_mask][0]

    eqw, eqwerr, col, colerr, flag, velcent, velwidth = eqwrange(ion, wl, flux, ferr,
                                                                 vrange, w_obs, f0, silent=silent, \
                                                                 sat_limit = sat_limit, sig_limit = sig_limit,
                                                                 plots = plots, plot_dir = plot_dir)
    eqw_list.append(eqw)
    eqwerr_list.append(eqwerr)
    col_list.append(col)
    colerr_list.append(colerr)
    flag_list.append(flag)
    velcent_list.append(velcent)
    velwidth_list.append(velwidth)

    return eqw_list, eqwerr_list, col_list, colerr_list, flag_list, velcent_list, velwidth_list



def json_eqw(json_file, fits_file, outfile, overwrite = False):
        ion_name = []; restwave = []; col = []; col_err = []; flag = []; zcomp = [];
        if not os.path.isfile(json_file):
                return np.array(ion_name), np.array(restwave), np.array(col), \
                        np.array(col_err), np.array(flag), np.array(zcomp)

        if not os.path.isfile(outfile) or overwrite == True:
                out =  open(outfile, 'w')
                out.write('#ion_name, restwave, logN, logN_err, flag_N, zcomp\n')
                
                with open(json_file) as data_file:
                        igmg_dict = json.load(data_file)
                        
                comp_list = []
                for ii, key in enumerate(igmg_dict['cmps'].keys()):
                        comp = AbsComponent.from_dict(igmg_dict['cmps'][key],
                                                      chk_sep=False, chk_data=False, chk_vel=False)
                        comp_list += [comp]
                        
                spec = readspec(fits_file)
                for i, cc in enumerate(comp_list):
                        for j, al in enumerate(cc._abslines):
                                al.analy['spec']=spec
                                al.measure_ew()
                                al.measure_aodm()
                                
                                ion_name.append(al.ion_name)
                                restwave.append(al.wrest.value)
                                col.append(al.attrib['logN'])
                                sig_logN = al.attrib['sig_logN']
                                if np.isnan(sig_logN):
                                        sig_logN = 0.0
                                col_err.append(sig_logN)
                                flag.append(al.attrib['flag_N'])
                                zcomp.append(al.z)
                                out.write('%s\t%f\t%e\t%e\t%d\t%0.6f\n'%(al.ion_name, al.wrest.value, \
                                        al.attrib['logN'], sig_logN, al.attrib['flag_N'], al.z))
                out.close()
        if os.path.isfile(outfile) and len(open(outfile).readlines()) > 1:
                ion_name = np.loadtxt(outfile, unpack = True, skiprows = 1, usecols = 0, dtype = 'str')
                restwave, col, col_err, flag, zcomp = \
                        np.loadtxt(outfile, unpack = True, skiprows = 1, usecols =  (1, 2, 3, 4, 5))

        return np.array(ion_name), np.array(restwave), np.array(col), np.array(col_err), np.array(flag), np.array(zcomp)


def load_veeper_fit(veeper_fn):
        if os.path.isfile(veeper_fn):
                restwaves, cols, sigcols, bvals, sigbvals, vels, sigvels, zcomps = \
                        np.loadtxt(veeper_fn, unpack=True, skiprows = 1,
                                   usecols = (1,3,4,5,6,7,8, 18), delimiter = '|')
                veeper_ions, labels = np.loadtxt(veeper_fn, unpack=True, skiprows = 1, \
                                                 usecols = (19, 21), dtype = 'str', delimiter = '|')
                
                if type(veeper_ions) == np.str_:
                        veeper_ions = np.array([veeper_ions]); restwaves =  np.array([restwaves]); 
                        cols =  np.array([cols]); sigcols =  np.array([sigcols]); 
                        bvals =  np.array([bvals]); sigbvals =  np.array([sigbvals]);
                        vels =  np.array([vels]); sigvels =  np.array([sigvels]); 
                        labels =  np.array([labels]); zcomps = np.array([zcomps]);
                        
                for i in range(len(veeper_ions)):
                        temp    = veeper_ions[i]
                        veeper_ions[i] = temp.replace(" ", "")
                        
                if sigcols.size > 1:
                        mask = sigcols > 0
                        veeper_ions = veeper_ions[mask]
                        restwaves = restwaves[mask]
                        cols = cols[mask]
                        sigcols = sigcols[mask]
                        bvals = bvals[mask]
                        sigbvals = sigbvals[mask]
                        vels = vels[mask]
                        sigvels = sigvels[mask]
                        labels = labels[mask]
                        zcomps = zcomps[mask]
        else:
                restwaves = []; cols = []; sigcols = []; bvals = []; sigbvals = [];
                vels      = []; sigvels = []; flags = []; veeper_ions = []; labels = []; zcomps = [];
        return veeper_ions, restwaves, cols, sigcols, bvals, sigbvals, vels, sigvels, labels, zcomps
