import numpy as np
import emcee
from astropy.modeling.models import Voigt1D
import corner
import matplotlib.pylab as plt
import seaborn as sns
import spec_helper_functions as shf

def deriv(x):
    return np.diff(np.append(x, x[-1]))

def guess_parameters(vel, flux, fluxerr, sig_lim = 2):

    #x0 = np.median(vel[flux < 1.0 - 2*std])
    
    amp_guess = 1.0 - np.min(flux)
    
    dwave = deriv(vel)
    d_eqw = dwave * (1. - flux) 
    eqw = np.sum(d_eqw)
    eqw_err = np.sqrt(np.sum(deriv(vel)**2  * fluxerr)**2 )

    f_eqw = np.cumsum(d_eqw) / np.max(np.cumsum(d_eqw)) #cumulative eqw fraction                                   
    # interpolate for velocity at which fraction of eqw = 0.5                                                      
    # we'll take this as the velocity centroid                                                                     
    # trying something new:                                                                                        
    vcent = np.interp(0.5, f_eqw, vel)

    lnf_guess = np.log(0.5)

    # if it's a non-detection, set the amplitude guess to 0
    if eqw < sig_lim*eqw_err:
        amp_guess = 0


    fwhl_guess = eqw
    fwhg_guess = eqw
    sat_flux = flux[flux < 0.1]
    if len(sat_flux) > 3:
        amp_guess = 2 + eqw / 50.
        fwhl_guess /= 4
        fwhg_guess = eqw - fwhl_guess
    
    return [vcent, amp_guess, fwhl_guess, fwhg_guess, lnf_guess]

def normalize_flux(ion, wl, flux, ferr, redshift = 0, vmin = -150, vmax = 150):
    w0 = shf.restwave(ion, redshift)
    vel = (wl-w0)/w0*3e5
    qq = ((vel > -600) & (vel < vmin) | ((vel > vmax) & (vel < 600)))
    med = np.median(flux[qq])
    flux /= med
    ferr /= med
    return flux, ferr
                        
def lnlike(theta, x, y, yerr, sat_lim = 0.1):
    x0, amp, fwhml, fwhmg, lnf = theta
    model = voigt_fit(x, theta, sat_lim = sat_lim)
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprior(theta):
    x0, amp, fwhml, fwhmg, lnf = theta
    if (-250 < x0 < 250) and (0.0 <= amp < 10.0) and (0.0 < fwhml < 250.0)\
            and (0.0 < fwhmg < 250.0) and (-10.0 < lnf <= 1.0):
        return 0.0
    return -np.inf

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)


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


def fit_single_voigt_profile(x, y, yerr, initial_guess, nwalkers = 100, \
                                 niterations = 100, jump_size = 1e-2, corner_plot = False):
    ndim = 5
    pos = [initial_guess + jump_size*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
    sampler.run_mcmc(pos, niterations)
    samples = sampler.chain[:, int(nwalkers/2):, :].reshape((-1, ndim))
    samples[:, 4] = np.exp(samples[:, 4])
    w0_mcmc, amp_mcmc, fwhml_mcmc, fwhmg_mcmc, lnf_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),\
                            zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    
    if corner_plot:
        fig = corner.corner(samples, labels=["$w0$", "$amp$", "fwhm_L", "fwhm_g", "$\ln\,f$"],
                      truths=initial_guess)
        plt.savefig('mcmc_voigt_fit.png')
    
    return w0_mcmc, amp_mcmc, fwhml_mcmc, fwhmg_mcmc, lnf_mcmc


def plot_ion_fit(ion, x, y, yerr, theta, initial_guess = None, color = 'red', label = None, ax = None, sat_lim = 0.1, linestyle = 'solid'):
    if ax == None:
        fig, ax = plt.subplots(nrows=1, ncols=1)

    if initial_guess: 
        model = voigt_fit(x, initial_guess, sat_lim = sat_lim)
        ax.plot(x, model, label = 'first guess', color = 'black', linestyle = 'dashed', alpha = 0.6)

    model = voigt_fit(x, theta, sat_lim = sat_lim)
    print(ion, theta)
    ax.plot(x,y, color = 'black', alpha = 0.2)
    ax.errorbar(x,y,yerr, color = 'black', alpha = 0.1)
    ax.plot(x, model, color = color, label = label, linestyle = linestyle, linewidth = 2.5, zorder= 10)
    ax.annotate(ion, xy=(-150, 0), fontsize = 12)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlim(min(x), max(x))
#    ax.legend()

def plot_veeper_ion_fit(ax, ion_list, model, orientation, radius, time = 11.2, vmin = -300, \
                            vmax = 300, work_dir = '../../data/analyzed_spectra'):
    if time == 11.2:
        redshift = 0.2
    else:
        redshift = 0

    base = shf.spec_base_filename(orientation, model, time, radius)
    fit = '%s/%s/FitInspection.fits'%(work_dir,base)

    wl, flux, ferr = shf.load_spec_from_fits(fit)

    nrows = len(ax)
    ncols = len(ax[0])

    for i, ion in enumerate(ion_list):
        row = int(i / ncols)
        col = i - row*ncols
        vv_ion, flux_ion, ferr_ion = ion_velocity_range(ion, wl, flux, ferr, redshift = redshift, vmin=vmin, vmax = vmax)
        ax[row][col].plot(vv_ion, flux_ion, color = 'green',  alpha = 0.7, linewidth = 2.5, label = 'veeper')
    


def voigt_fit(x, theta, sat_lim):
    model = 1.0 - Voigt1D(x_0=theta[0], amplitude_L=theta[1], fwhm_L=theta[2], fwhm_G=theta[3])(x)
    sat= model< sat_lim
    model[sat] = sat_lim
    return model
    

# takes in ion wavelength, flux, flux error and 
# returns the velocity range +- 200 km/s around the central velocity
def ion_velocity_range(ion, wl, flux, fluxerr, redshift = 0, vmin = -200, vmax = 200, ):
    w0 = shf.restwave(ion, redshift) # REST wavelength for rest-frame EW; observed wavelength for observed EW 
    vv = (wl-w0) / w0 * 2.9979e5
    mask = (vv > vmin) & (vv < vmax)
    return vv[mask], flux[mask], fluxerr[mask]


 
def fit_spectrum(model, orientation, radius, ion_list = [], time = 11.2, normalize = False,\
                vmin=-150, vmax =150, nwalkers = 100, niterations = 100, save_fit = True, \
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
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*nrows, 3.8*ncols))
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
        initial_guess = guess_parameters(vv_ion, flux_ion, ferr_ion)
        result = fit_single_voigt_profile(vv_ion, flux_ion, ferr_ion, \
                                          initial_guess, niterations=niterations);
        theta = [result[0][0], result[1][0], result[2][0], result[3][0]]

        plot_ion_fit(ion, vv_ion, flux_ion, ferr_ion, theta, initial_guess = initial_guess, \
                         color = 'red', ax = ax, label = "mcmc fit" )

        vv_list.append(vv_ion)
        flux_list.append(flux_ion)
        ferr_list.append(ferr_ion)
        theta_list.append(ferr_ion)
        if save_fit:
            outfile.write('%s %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n'\
                              %(ion, result[0][0], result[0][1], result[0][2], \
                                    result[1][0], result[1][1], result[1][2],\
                                    result[2][0], result[2][1],result[2][2], \
                                    result[3][0], result[3][1],result[3][2], \
                                    result[4][0], result[4][1],result[4][2]))
#    plt.show()
    if plot_veeper:
        plot_veeper_ion_fit(axes, ion_list, model, orientation, radius, vmin = vmin, vmax = vmax)
#    ax[0][0].legend()#loc = 'lower right')
    plt.savefig('../../plots/mcmc/%s_%s_%ikpc.png'%(orientation, model, radius))
    return vv_list, flux_list, ferr_list, theta_list
