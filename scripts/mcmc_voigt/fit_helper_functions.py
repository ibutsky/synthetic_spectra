import numpy as np
import emcee
import sys
from astropy.modeling.models import Voigt1D
import corner
import matplotlib.pylab as plt
import seaborn as sns
import spec_helper_functions as shf
from scipy import signal


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

def guess_parameters(vel, flux, fluxerr, sig_lim = 2, sat_lim = 0.05):

    amp_guess = 1.0 - np.min(flux)
    if amp_guess < 0:
        amp_guess = 0

    eqw, eqw_err = calc_eqw(vel, flux, fluxerr)
    d_eqw = deriv(vel) * (1. - flux)
    f_eqw = np.cumsum(d_eqw) / np.max(np.cumsum(d_eqw)) 
    lnf_guess = np.log(0.5)

    # if it's a non-detection, set the amplitude guess to 0
    if eqw < sig_lim*eqw_err:
        amp_guess = 0

    fwhl_guess = eqw 
    fwhg_guess = eqw 
    sat_flux = flux[flux < sat_lim]
    if len(sat_flux) > 3:
        amp_guess =  eqw / 40.
        fwhl_guess = 50 + eqw / 50;
        fwhg_guess = fwhl_guess
    
    vcent_peaks = find_num_peaks(vel, flux)
    if len(vcent_peaks) == 0:
        vcent = 0
        amp_guess = 0
    else:
        vcent = vcent_peaks[0]
    theta = [vcent, amp_guess, fwhl_guess, fwhg_guess, lnf_guess]
    if len(vcent_peaks) > 1:
        theta = [vcent_peaks[0], amp_guess, fwhl_guess/2., fwhg_guess/2., \
                 vcent_peaks[1], amp_guess, fwhl_guess/2., fwhg_guess/2.,lnf_guess]


    return theta

def normalize_flux(ion, wl, flux, ferr, redshift = 0, vmin = -150, vmax = 150):
    w0 = shf.restwave(ion, redshift)
    vel = (wl-w0)/w0*3e5
    qq = ((vel > -600) & (vel < vmin) | ((vel > vmax) & (vel < 600)))
    med = np.median(flux[qq])
    flux /= med
    ferr /= med
    return flux, ferr
                        
def lnlike(theta, x, y, yerr, sat_lim = 0.05):
#    x0, amp, fwhml, fwhmg, x02, amp2, fwhml2, fwhmg2, lnf = theta
    lnf = theta[-1]
    model = voigt_fit(theta, x, y, yerr, sat_lim = sat_lim)
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def calc_prior(x0, amp, fwhml, fwhmg, lnf, vmin = -250, vmax = 250):
    prior = -np.inf
    if (vmin < x0 < vmax) and (0.0 <= amp    < 10.00) and (0.0 <= fwhml  < 250.0) \
                          and (0.0 <= fwhmg  < 250.0) and (-10. < lnf    <=  1.0):
        prior = 0.0
    return prior 

def lnprior(theta, x, y, yerr, sig_lim = 2, sat_lim = 0.05,  vmin = -250, vmax = 250):
    #todo don't need x, y, yerr or siglim in this
    prior = calc_prior(theta[0], theta[1], theta[2], theta[3], theta[-1])
    if len(theta) == 9:
        prior += calc_prior(theta[4], theta[5], theta[6], theta[7], theta[-1])
    return prior

def lnprob(theta, x, y, yerr, sig_lim = 2, sat_lim = 0.05):
    lp = lnprior(theta, x, y, yerr, sig_lim = sig_lim, sat_lim = sat_lim)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr, sat_lim = sat_lim)


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


def fit_single_voigt_profile(x, y, yerr, initial_guess, nwalkers = 100,\
                                 niterations = 100, jump_size = 1e-2, corner_plot = False):

    ndim  = len(initial_guess)
    pos = [initial_guess + jump_size*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
    sampler.run_mcmc(pos, niterations)
    samples = sampler.chain[:, int(nwalkers/2):, :].reshape((-1, ndim))
    samples[:, ndim-1] = np.exp(samples[:, ndim-1])
    if ndim == 5:
        w0_mcmc, amp_mcmc, fwhml_mcmc, fwhmg_mcmc, lnf_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),\
                            zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        result = [w0_mcmc, amp_mcmc, fwhml_mcmc, fwhmg_mcmc, lnf_mcmc]
        theta = [result[0][0], result[1][0], result[2][0], result[3][0], result[4][0]]
    elif ndim == 9:
        w0_mcmc, amp_mcmc, fwhml_mcmc, fwhmg_mcmc, w02_mcmc, amp2_mcmc, \
            fwhml2_mcmc, fwhmg2_mcmc, lnf_mcmc  = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),\
                            zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        result = [w0_mcmc, amp_mcmc, fwhml_mcmc, fwhmg_mcmc, w02_mcmc, amp2_mcmc, fwhml2_mcmc, fwhmg2_mcmc, lnf_mcmc]
        theta = [result[0][0], result[1][0], result[2][0], \
                 result[3][0], result[4][0], result[5][0], \
                 result[6][0], result[7][0], result[8][0]]
    
    if corner_plot:
        if len(initial_guess) == 5:
            labels = ["$w0$", "$amp$", "fwhm_L", "fwhm_g", "$\ln\,f$"]
        elif len(initial_guess) == 9:
            labels = ["$w0$", "$amp$", "fwhm_L", "fwhm_g", "$w02$", "$amp2$", "fwhm_L2", "fwhm_g2","$\ln\,f$"]
        fig = corner.corner(samples, labels=labels, truths=theta)
        plt.savefig('mcmc_voigt_fit.png')

    return result


def plot_ion_fit(ion, x, y, yerr, theta, initial_guess = None, color = 'red', label = None, ax = None, sat_lim = 0.05, linestyle = 'solid'):
    if ax == None:
        fig, ax = plt.subplots(nrows=1, ncols=1)

    if initial_guess: 
        model = voigt_fit(initial_guess, x, y, yerr, sat_lim = sat_lim)
        ax.plot(x, model, label = 'first guess', color = 'black', linestyle = 'dashed', alpha = 0.6)

    model = voigt_fit(theta, x, y, yerr, sat_lim = sat_lim)
    ax.plot(x,y, color = 'black', alpha = 0.2)
    ax.errorbar(x,y,yerr, color = 'black', alpha = 0.1)
    ax.plot(x, model, color = color, label = label, linestyle = linestyle, linewidth = 2.5, zorder= 10)
    ax.annotate(ion, xy=(-150, 0), fontsize = 12)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlim(min(x), max(x))
#    ax.legend()

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

    vv_ion, flux_ion, ferr_ion = ion_velocity_range(ion, wl, flux, ferr, redshift = redshift, vmin=vmin, vmax = vmax)
    ax.plot(vv_ion, flux_ion, color = 'green',  alpha = 0.7, linewidth = 2.5, label = 'veeper')
    


def voigt_fit(theta, x, y, yerr, sat_lim = 0.05):
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
    return vv[mask], flux[mask], fluxerr[mask]


 
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
        vv_ion, flux_ion, ferr_ion = ion_velocity_range(ion, wl, flux, ferr, redshift = redshift, vmin=vmin, vmax = vmax)
        initial_guess = guess_parameters(vv_ion, flux_ion, ferr_ion, sat_lim = sat_lim, sig_lim = sig_lim)

        result = fit_single_voigt_profile(vv_ion, flux_ion, ferr_ion,  initial_guess, \
                                                     corner_plot = corner_plot, niterations=niterations);
       
        if len(result) == 5:
            theta = [result[0][0], result[1][0], result[2][0], \
                 result[3][0], result[4][0]]

        elif len(result) == 9:
            theta = [result[0][0], result[1][0], result[2][0], \
                     result[3][0], result[4][0], result[5][0], \
                     result[6][0], result[7][0], result[8][0]]    

#        prob_single = lnprob(theta_single, vv_ion, flux_ion, ferr_ion)
#        prob_double = lnprob(theta_double, vv_ion, flux_ion, ferr_ion)
#        dL_single = result_single[2][2] - result_single[2][1]
#        dL_double = result_double[2][2]- result_double[2][1]
#        print(ion, prob_single, prob_double, theta_single[-1], theta_double[-1])
#        print(ion, result_single[2][0], dL_single, result_double[2][0], dL_double)
#        print(theta_single, theta_double)
#        if prob_single > prob_double:
#            theta = theta_single
#            initial_guess = initial_guess_single
#            result = result_single
#        else:
#            theta = theta_double
#            initial_guess = initial_guess_double
#            result = result_double

        if plot_veeper:
            plot_veeper_ion_fit(ion, model, orientation, radius, ax = ax, vmin = vmin, vmax = vmax)
        plot_ion_fit(ion, vv_ion, flux_ion, ferr_ion, theta, initial_guess = initial_guess, sat_lim = sat_lim, \
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
                                    result[3][0], result[3][1],result[3][2], \
                                    result[4][0], result[4][1],result[4][2]))
        if row == nrows-1:
            ax.set_xlabel('Velocity (km/s)')
        if col == 0:
            ax.set_ylabel('Flux')
        if i == 0:
            ax.legend(frameon = True)
           
    plt.savefig('../../plots/mcmc/%s_%s_%ikpc.png'%(orientation, model, radius))
    return vv_list, flux_list, ferr_list, theta_list
