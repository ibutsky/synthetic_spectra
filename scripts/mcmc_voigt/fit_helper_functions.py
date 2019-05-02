import numpy as np
import emcee
import sys
from astropy.modeling.models import Voigt1D
import corner
import matplotlib.pylab as plt
import seaborn as sns
import spec_helper_functions as shf
from scipy import signal

from scipy.special import wofz
from astropy.constants import c as C
echarge = 4.803204505713468e-10
m_e = 9.10938291e-28
c = C.to('cm/s').value
c2 = 29979245800.0

def Hfunc(x,a):
    z=x+1j*a
    I=wofz(z).real
    return I

def find_ion_props(ion, redshift = 0):
    ion = ion.replace(" ", "")
    if ion == 'HI':
        restwave = 1215.67
        gamma =  4.690000e+08 
        fosc =  4.160000e-01
    elif ion == 'CII':
        restwave = 1334.5323
        gamma = 2.410000e+08    
        fosc = 1.290000e-01
    elif ion == 'CIII':
        restwave = 977.020
        gamma = 1.767000e+09
        fosc = 7.586e-1
    elif ion == 'CIV':
        restwave = 1548.204
        gamma = 2.650000e+08    
        fosc = 1.900000e-01
    elif ion == 'NV':
        restwave = 1238.821
        gamma = 3.400000e+08    
        fosc = 1.560000e-01
    elif ion == 'NIII':
        restwave = 989.799
        gamma = 4.180000e+08 
        fosc = 1.230000e-01
    elif ion == 'SiII':
        restwave = 1190.4158
        gamma =  6.530000e+08
        fosc = 2.770000e-01
    elif ion == 'SiIII':
        restwave = 1206.5
        gamma = 2.480000e+09    
        fosc = 1.630000e+00
    elif ion == 'SiIV':
        restwave = 1402.7729
        gamma = 8.630000e+08    
        fosc = 2.550000e-01
    elif ion == 'MgII':
        restwave = 1239.9253
    elif ion == 'OVI':
        restwave = 1031.9261
        gamma = 4.160000e+08    
        fosc = 1.330000e-01
    else:
        print("%s not recognized as ion"%(ion))
    return restwave, gamma, fosc


def veeper_voigt_fit(waves, theta, ion_props, z = 0, sat_lim = 0.1):
    peak1 = 1-single_voigt_fit(waves, theta[0], theta[1], theta[2], ion_props)
    model = 1 - peak1
    if len(theta) > 4:
        peak2 = 1 - single_voigt_fit(waves, theta[3], theta[4], theta[5], ion_props)
        model = 1 - (peak1 + peak2)

    return model

def single_voigt_fit(waves, coldens, bval, vels, ion_props, z = 0, sat_lim = 0.1):
    tautot=np.zeros(len(waves))
    thatfactor=(1.+z)

    lam = waves
    lam0, gam, fosc = ion_props
    dlam=bval*lam0/2.9979e5  #Doppler param in wavelength                                                                                
    x=(lam-lam0-lam0*vels/2.9979e5)/dlam
    a=gam/(4.*np.pi*(c*1e13/lam0**2*dlam))
    vp=Hfunc(x,a)
    tauval= np.sqrt(np.pi) * echarge ** 2 / m_e / c ** 2 * (lam0 * 1e-8) ** 2 / dlam * 1e8 * (10 ** coldens) * fosc * vp
    tautot+=tauval
    model =  np.exp(-tautot)

    sat = model < sat_lim
    model[sat] = np.mean(sat_lim)

    return model

def deriv(x):
    return np.diff(np.append(x, x[-1]))

def smooth(x, window_len = 30, window = 'blackman'):
    if x.ndim != 1:
        print("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        print("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x
       
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        print("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')     
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len-1:-window_len+1] 


def find_num_peaks(vv_ion, flux_ion, vcent = 0, width = 100, window_len = 30, window = 'blackman'):    
    # Finds the number of peaks in the absorption spectrum to decide wether to fit one line or two
    fsmooth = smooth(flux_ion, window_len = window_len, window = window)
    vsmooth = smooth(vv_ion, window_len = window_len, window = window)
    peaks = signal.find_peaks(1/fsmooth, height = 1.1)[0]

    return vsmooth[peaks]

def calc_eqw(vel, flux, fluxerr):
    dwave = deriv(vel)
    d_eqw = dwave * (1. - flux)
    eqw = np.sum(d_eqw)
    eqw_err = np.sqrt(np.sum((deriv(vel)**2  * fluxerr)**2) )

    return eqw, eqw_err

def guess_parameters(waves, flux, fluxerr, ion_props, sig_lim = 2, sat_lim = 0.1):
    w0 = np.median(waves)
    vel = (waves - w0) / w0*2.997e5
    amp_guess = 1.0 - np.min(flux)
    if amp_guess < 0:
        col_guess = 0
        amp_guess = 0
    col_guess = 12 + 2*amp_guess

    eqw, eqw_err = calc_eqw(vel, flux, fluxerr)
    d_eqw = deriv(vel) * (1. - flux)
    f_eqw = np.cumsum(d_eqw) / np.max(np.cumsum(d_eqw)) 
    lnf_guess = np.log(0.5)

    # if it's a non-detection, set the amplitude guess to 0
    if eqw < sig_lim*eqw_err:
        col_guess = 0

    
    bval_guess = eqw
    sat_flux = flux[flux < sat_lim]
    if len(sat_flux) > 3:
        col_guess =  16.
        bval_guess = 30 + eqw / 50;

    
    vcent_peaks = find_num_peaks(vel, flux)
    if len(vcent_peaks) == 0:
        vcent = 0
        col_guess = 0
    else:
        vcent = vcent_peaks[0]

    theta = [col_guess, bval_guess, vcent, lnf_guess]

    if len(vcent_peaks) > 1:
        vcent2 = vcent_peaks[1]
        bval_guess= np.abs(vcent2 - vcent)/2
        theta = [col_guess + 1, bval_guess, vcent,\
                 col_guess + 1, bval_guess, vcent2, lnf_guess]
    print(theta)    
    return theta

def normalize_flux(ion, wl, flux, ferr, redshift = 0, vmin = -150, vmax = 150):
    w0 = shf.restwave(ion, redshift)
    vel = (wl-w0)/w0*3e5
    qq = ((vel > -600) & (vel < vmin) | ((vel > vmax) & (vel < 600)))
    med = np.median(flux[qq])
    flux /= med
    ferr /= med
    return flux, ferr
                        
def lnlike(theta, x, y, yerr, ion_props, sat_lim = 0.1):
    lnf = theta[-1]
    model = veeper_voigt_fit(x, theta, ion_props)
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return  -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def calc_prior(col, bval, x0, lnf, vmin = -250, vmax = 250):
    prior = -np.inf
    if (vmin < x0 < vmax) and (0 <  col < 22) and (0.0 <  bval < 250.0) and (-10. < lnf    <=  1.0):
        prior = 0.0
    return prior 

def lnprior(theta):
    #todo don't need x, y, yerr or siglim in this
    prior = calc_prior(theta[0], theta[1], theta[2], theta[-1])
    if len(theta) == 7:
        prior += calc_prior(theta[3], theta[4], theta[5], theta[-1])
    return prior

def lnprob(theta, x, y, yerr, ion_props, sig_lim = 2, sat_lim = 0.1):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr, ion_props, sat_lim = sat_lim)


def set_rows_cols(num_plots):
    if num_plots < 1:
        nrows = 0
        ncols = 0
    if num_plots == 1:
        nrows = 1
        ncols = 1
    elif num_plots < 9:
        nrows = 2
        ncols = int(num_plots / 2.)
    else:
        nrows = 3
        ncols = int(num_plots / 3.)
    return nrows, ncols


def fit_single_voigt_profile(ion_props, x, y, yerr, initial_guess, nwalkers = 100,\
                                 niterations = 100, jump_size = 1e-3, corner_plot = False):

    ndim  = len(initial_guess)
    pos = [initial_guess + jump_size*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr, ion_props))
    sampler.run_mcmc(pos, niterations)
    discard = int(nwalkers/2)
    discard = 499
    samples = sampler.chain[:, :, :].reshape((-1, ndim))
    samples[:, ndim-1] = np.exp(samples[:, ndim-1])
    if ndim == 4:
        col, bval, vel, lnf = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),\
                            zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        result = [col, bval, vel, lnf]
        theta = [result[0][0], result[1][0], result[2][0], result[3][0]]
    elif ndim == 7:
        col, bval, vel, col2, bval2, vel2, lnf  = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),\
                            zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        result = [col, bval, vel, col2, bval2, vel2, lnf]
        theta = [result[0][0], result[1][0], result[2][0], \
                 result[3][0], result[4][0], result[5][0], \
                 result[6][0]]
    
    if corner_plot:
        if len(initial_guess) == 4:
            labels = ["N", "bval","vel", "ln"]
        if len(initial_guess) == 7:
            labels = ["N", "bval","vel","N2", "bval2", "vel2",  "ln"]
        fig = corner.corner(samples, labels=labels, truths=theta)
        plt.savefig('mcmc_voigt_fit.png')
    return result


def plot_ion_fit(ion, wl, vel,  y, yerr, theta, ion_props, initial_guess = None, color = 'red', label = None, ax = None, sat_lim = 0.05, linestyle = 'solid'):
    if ax == None:
        fig, ax = plt.subplots(nrows=1, ncols=1)
    if initial_guess: 
        model = veeper_voigt_fit(wl, initial_guess, ion_props, z = 0.2)
        ax.plot(vel, model, label = 'first guess', color = 'black', linestyle = 'dashed', alpha = 0.6)

    model = veeper_voigt_fit(wl, theta,ion_props,  z = 0.2)
    ax.plot(vel,y, color = 'black', alpha = 0.2)
    ax.errorbar(vel,y,yerr, color = 'black', alpha = 0.1)
    ax.plot(vel, model, color = color, label = label, linestyle = linestyle, linewidth = 2.5, zorder= 10)
    ax.annotate(ion, xy=(-150, 0), fontsize = 12)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlim(min(vel), max(vel))

def plot_veeper_ion_fit(ion, model, orientation, radius, ax=None, time = 11.2, vmin = -300, \
                            vmax = 300, work_dir = '../../data/analyzed_spectra'):

    if ax == None:
        fig, ax = plt.subplots(nrows=1, ncols=1)
    
    if time == 11.2:
        redshift = 0.2
    else:
        redshift = 0

    base = shf.spec_base_filename(orientation, model, time, radius)
    fit = '%s/%s/FitInspection.fits'%(work_dir,base)

    wl, flux, ferr = shf.load_spec_from_fits(fit)

    wl_ion, vv_ion, flux_ion, ferr_ion = ion_velocity_range(ion, wl, flux, ferr, redshift = redshift, vmin=vmin, vmax = vmax)

    ax.plot(vv_ion, flux_ion, color = 'green',  alpha = 0.7, linewidth = 2.5, label = 'veeper')
    


def voigt_fit(theta, x, y, yerr, sat_lim = 0.1):
    peak1 = Voigt1D(x_0=theta[0], amplitude_L=theta[1], fwhm_L=theta[2], fwhm_G=theta[3])(x)
    peak2 = np.zeros(len(peak1))
    if len(theta) > 5:
        peak2 = Voigt1D(x_0=theta[4], amplitude_L=theta[5], fwhm_L=theta[6], fwhm_G=theta[7])(x)
    model = 1.0 - (peak1 + peak2)
    sat = model < sat_lim
    model[sat] = np.mean(sat_lim)
    return model
    

# takes in ion wavelength, flux, flux error and 
# returns the velocity range +- 200 km/s around the central velocity
def ion_velocity_range(ion, wl, flux, fluxerr, redshift = 0, vmin = -200, vmax = 200, ):
    w0 = shf.restwave(ion, redshift) # REST wavelength for rest-frame EW; observed wavelength for observed EW 
    vv = (wl-w0) / w0 * 2.9979e5
    mask = (vv > vmin) & (vv < vmax)
    wl = wl / (1+ redshift)
    return wl[mask], vv[mask], flux[mask], fluxerr[mask]


 
def fit_spectrum(model, orientation, radius, ion_list = [], time = 11.2, normalize = False, corner_plot = False, \
                 vmin=-150, vmax =150, nwalkers = 100, niterations = 100, save_fit = True, sat_lim = 0.05, sig_lim = 0.05, \
                 work_dir = '../../data/analyzed_spectra', use_errors = True, plot_veeper = True):

    if time == 11.2:
        redshift = 0.2
    else:
        redshift = 0
        print('ARE YOU SURE YOU WANT REDSHIFT = 0?')

    base = shf.spec_base_filename(orientation, model, time, radius)
    wl, flux, ferr = shf.load_spec_from_fits('%s/%s/%s_ibnorm.fits'%(work_dir,base, base))
    if not use_errors:
        ferr = np.zeros(len(flux))

    if save_fit:
        outfile = open('%s/%s/mcmc_fit.dat'%(work_dir, base), 'w')
        outfile.write('ion_name v_center[1] verr1 verr2 amp[4] amperr1 amperr2 fwhm_l[7] lerr1 \
lerr2 fwhm_g[10] gerr1 gerr2 lnf[13] lnferr1 lnferr2\n')
    # plotting is temporary
    ncols = int(len(ion_list) / 2)
    nrows, ncols = set_rows_cols(len(ion_list))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 3.8*nrows))
    vv_list = []
    flux_list = []
    ferr_list = []
    theta_list = []
    for i, ion in enumerate(ion_list):
        row = int(i / ncols)
        col = i - row*ncols
        if nrows == 1:
            if ncols == 1:
                ax = axes
            else:
                ax = axes[i]
        else:
            ax = axes[row][col]

        if normalize:
            flux, ferr =  normalize_flux(ion, wl, flux, ferr, redshift = redshift, vmin = vmin, vmax = vmax)
        wl_ion, vv_ion, flux_ion, ferr_ion = ion_velocity_range(ion, wl, flux, ferr, redshift = redshift, vmin=vmin, vmax = vmax)
        ion_props = find_ion_props(ion, redshift = redshift)
        initial_guess = guess_parameters(wl_ion, flux_ion, ferr_ion, ion_props, sat_lim = sat_lim, sig_lim = sig_lim)

        result = fit_single_voigt_profile(ion_props, wl_ion, flux_ion, ferr_ion,  initial_guess, \
                                                     corner_plot = corner_plot, niterations=niterations);
        if len(result) == 4:
            theta = [result[0][0], result[1][0], result[2][0], result[3][0]]

        elif len(result) == 7:
            theta = [result[0][0], result[1][0], result[2][0], \
                     result[3][0], result[4][0], result[5][0], \
                     result[6][0]]    


        if plot_veeper:
            plot_veeper_ion_fit(ion, model, orientation, radius, ax = ax, vmin = vmin, vmax = vmax)
        plot_ion_fit(ion, wl_ion, vv_ion, flux_ion, ferr_ion, theta, ion_props, initial_guess = initial_guess, sat_lim = sat_lim, \
                         color = 'red', ax = ax, label = "mcmc fit" )

        vv_list.append(vv_ion)
        flux_list.append(flux_ion)
        ferr_list.append(ferr_ion)
        theta_list.append(theta)
        if save_fit:
            outfile.write('%s %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n'\
                              %(ion, result[0][0], result[0][1], result[0][2], \
                                    result[1][0], result[1][1], result[1][2],\
                                    result[2][0], result[2][1],result[2][2], \
                                    result[3][0], result[3][1],result[3][2]))#, \
                           #         result[4][0], result[4][1],result[4][2]))
        if row == nrows-1:
            ax.set_xlabel('Velocity (km/s)')
        if col == 0:
            ax.set_ylabel('Flux')
        if i == 0:
            ax.legend(frameon = True)
           
#    plt.savefig('../../plots/mcmc/%s_%s_%ikpc.png'%(orientation, model, radius))
    return vv_list, flux_list, ferr_list, theta_list
