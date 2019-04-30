import numpy as np
import emcee
from astropy.modeling.models import Voigt1D
import corner
import matplotlib.pylab as plt
import seaborn as sns
import spec_helper_functions as shf


def guess_parameters(x, flux):
    min_x = np.min(x)
    min_flux = np.min(flux)
    inv_flux = 1.0 - flux
    max_flux = np.max(inv_flux)
    center_x_range = x[inv_flux > 0.7 * inv_flux]
    x0 = np.median(center_x_range)
    
    #x0 = x[flux == min_flux][0]
    amp_guess = 1.0 + min_flux
    hm =0.5*(1.0 + min_flux)

    fwhm_guess = np.interp(hm, flux[x > x0], x[x > x0]) - x0
    lnf_guess = np.log(0.5)
    if fwhm_guess < 10:
        amp_guess = 0
    return [x0, amp_guess, fwhm_guess, fwhm_guess, lnf_guess]

def normalize_flux(ion, wl, flux, ferr,vmin = -150, vmax = 150):
    w0 = shf.restwave(ion)
    vel = (wl-w0)/w0*3e5
    qq = ((vel > -600) & (vel < vmin) | ((vel > vmax) & (vel < 600)))
    med = np.median(flux[qq])
    flux /= med
    ferr /= med
    return flux, ferr
                        
def lnlike(theta, x, y, yerr):
    x0, amp, fwhml, fwhmg, lnf = theta
    v1 = Voigt1D(x_0 = x0, amplitude_L = amp, fwhm_L = fwhml, fwhm_G = fwhmg)
    model = 1.0 - v1(x)
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprior(theta):
    x0, amp, fwhml, fwhmg, lnf = theta
    
    if (-250 < x0 < 250) and (0.0 <= amp < 5.0) and (0.0 < fwhml < 250.0)\
            and (0.0 < fwhmg < 250.0) and (-10.0 < lnf <= 1.0):
        return 0.0
    return -np.inf

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)


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




def plot_ion_fit(ion, x, y, yerr, theta, initial_guess = None, color = 'red', label = None, ax = None, linestyle = 'solid'):
    if ax == None:
        fig, ax = plt.subplots(nrows=1, ncols=1)

    if initial_guess: 
        v = Voigt1D(x_0=initial_guess[0], amplitude_L=initial_guess[1], fwhm_L=initial_guess[2], fwhm_G=initial_guess[3])
        ax.plot(x, 1.0-v(x), label = 'first guess', color = 'black', linestyle = 'dashed', alpha = 0.6)
    v1 = Voigt1D(x_0=theta[0], amplitude_L=theta[1], fwhm_L=theta[2], fwhm_G=theta[3])
    ax.plot(x,y, color = 'black', alpha = 0.2)
    ax.errorbar(x,y,yerr, color = 'black', alpha = 0.1)
    ax.plot(x, 1.0-v1(x), color = color, label = label, linestyle = linestyle, linewidth = 2.5)
    ax.annotate(ion, xy=(200, 0), fontsize = 12)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlim(min(x), max(x))
#    ax.legend()

def plot_veeper_ion_fit(ax, ion_list, model, orientation, radius, vmin = -300, vmax = 300):
    fn = shf.spec_folder(orientation, model)+'%ikpc_VPmodel.fits'%(radius)
    wl, flux, ferr = shf.load_spec_from_fits(fn)
    nrows = 2
    ncols = 3
    for i, ion in enumerate(ion_list):
        row = int(i / ncols)
        col = i - row*ncols
        vv_ion, flux_ion, ferr_ion = ion_velocity_range(ion, wl, flux, ferr, vmin=vmin, vmax = vmax)
        ax[row][col].plot(vv_ion, flux_ion, color = 'red', linewidth = 2.5, label = 'veeper')
    


def voigt_fit(x, theta):
    return 1.0 - Voigt1D(x_0=theta[0], amplitude_L=theta[1], fwhm_L=theta[2], fwhm_G=theta[3])(x)
    

# takes in ion wavelength, flux, flux error and 
# returns the velocity range +- 200 km/s around the central velocity
def ion_velocity_range(ion, wl, flux, fluxerr, vmin = -200, vmax = 200):
    w0 = shf.restwave(ion) # REST wavelength for rest-frame EW; observed wavelength for observed EW 
    vv = (wl-w0) / w0 * 2.9979e5
    mask = (vv > vmin) & (vv < vmax)
    return vv[mask], flux[mask], fluxerr[mask]

#def fit_ion_velocity(ion, vv, flux, fluxerr, vmin=-200, vmax = 200):
 
def fit_spectrum(model, orientation, radius, ion_list = [], \
                vmin=-150, vmax =150, nwalkers = 100, niterations = 100, save_fit = True, \
                 work_dir = '../../data/analyzed_spectra'):
        

    fn = shf.spec_base_filename(orientation, model, radius)+'_ibnorm.fits'
    wl, flux, ferr = shf.load_spec_from_fits(fn)
    if save_fit:
        outfile = open(shf.spec_folder(orientation, model)+'ibfit_%ikpc.dat'%(radius), 'w')
        outfile.write('ion_name v_center[1] verr1 verr2 amp[4] amperr1 amperr2 fwhm_l[7] lerr1 \
lerr2 fwhm_g[10] gerr1 gerr2 lnf[13] lnferr1 lnferr2\n')
    # plotting is temporary
    ncols = int(len(ion_list) / 2)
    fig, ax = plt.subplots(nrows=2, ncols=ncols, figsize=(12, 6))
    vv_list = []
    flux_list = []
    ferr_list = []
    theta_list = []
    for i, ion in enumerate(ion_list):
        row = int(i / ncols)
        col = i - row*ncols
        flux, ferr =  normalize_flux(ion, wl, flux, ferr, vmin = vmin, vmax = vmax)
        vv_ion, flux_ion, ferr_ion = ion_velocity_range(ion, wl, flux, ferr, vmin=vmin, vmax = vmax)
        initial_guess = guess_parameters(vv_ion, flux_ion)
        result = fit_single_voigt_profile(vv_ion, flux_ion, ferr_ion, \
                                          initial_guess, niterations=niterations);
        theta = [result[0][0], result[1][0], result[2][0], result[3][0]]

        plot_ion_fit(ion, vv_ion, flux_ion, ferr_ion, theta, initial_guess = initial_guess, \
                         color = 'green', ax = ax[row][col], label = "mcmc fit" )
        if row == 0 and col == 0:
             ax[row][col].legend(loc = 3)
#        plot_veeper_ion_fit(ion, model, orientation, radius, ax = ax[row][col], vmin = vmin, vmax = vmax)

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
#    ax[0][0].legend(loc = 3)
    plot_veeper_ion_fit(ax, ion_list, model, orientation, radius, vmin = vmin, vmax = vmax)
    plt.savefig('%s_%s_%ikpc.png'%(orientation, model, radius))
    return vv_list, flux_list, ferr_list, theta_list
