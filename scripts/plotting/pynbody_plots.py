import matplotlib
matplotlib.use('Agg')
import sys
import pynbody
import pynbody.plot as pp
import pynbody.plot.sph as sph
import matplotlib.pylab as plt
import numpy as np
import h5py as h5

pynbody.config['number_of_threads'] = 2



def make_plots():
    output = int(sys.argv[1])
    plot_type = sys.argv[2]
    cr = int(sys.argv[3])

    if cr:
        s = pynbody.load('/nobackup/ibutsky/simulations/patient0_agncr/pioneer50h243.1536gst1bwK1BH.%06d'%(output))
    else:
        s = pynbody.load('/nobackup/ibutsky/simulations/patient0_nocr/pioneer50h243.1536gst1bwK1BH.%06d'%(output))
#        s = pynbody.load('/nobackup/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.%06d'%(output))

    s.physical_units()

    if plot_type == 'density':
        make_density_plots(s, output, cr)
    elif plot_type == 'rotation':
        make_rotation_curves(s, output, cr)
    elif plot_type == 'temperature':
        make_temperature_plots(s, output, cr)
    elif plot_type == 'all':
        make_density_plots(s, output, cr)
        make_rotation_curves(s, output, cr)
    elif plot_type == 'sfh':
         make_sfh(s, output, cr)

def make_density_plots(s, output, cr):
    pynbody.analysis.angmom.faceon(s.g)
    sph.image(s.g,qty="rho",units="g cm^-3",width=60,cmap="magma", vmin=5e-27, vmax=1e-23)
    if cr:
        plt.savefig('pioneer.%06d_density_faceon.png'%(output))
    else:
        plt.savefig('pioneer.%06d_density_faceon_nocr.png'%(output))
    plt.clf()

    pynbody.analysis.angmom.sideon(s.s)
    sph.image(s.g,qty="rho",units="g cm^-3",width=60,cmap="magma", vmin=5e-27, vmax=1e-23)
    if cr:
        plt.savefig('pioneer.%06d_density_sideon.png'%(output))
    else:
        plt.savefig('pioneer.%06d_density_sideon_nocr.png'%(output))
    plt.clf()

def make_temperature_plots(s, output, cr):
    pynbody.analysis.angmom.faceon(s.g)
    sph.image(s.g,qty="temp",units="K",width=100,cmap="afmhot", vmin=1e4, vmax=1e6)
    if cr:
        plt.savefig('pioneer.%06d_temperature_faceon.png'%(output))
    else:
        plt.savefig('pioneer.%06d_temperature_faceon_nocr.png'%(output))
    plt.clf()

    pynbody.analysis.angmom.sideon(s.s)
    sph.image(s.g,qty="temp",units="K",width=100,cmap="afmhot", vmin=1e4, vmax=1e6)
    if cr:
        plt.savefig('pioneer.%06d_temperature_sideon.png'%(output))
    else:
        plt.savefig('pioneer.%06d_temperature_sideon_nocr.png'%(output))
    plt.clf()

def make_sfh(s, output, cr):
    #s_agn = pynbody.load('/nobackup/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.003456')
    s_agn = pynbody.load('/nobackup/ibutsky/simulations/patient0_nocr/pioneer50h243.1536gst1bwK1BH.%06d'%(output))
#    s_sn = pynbody.load('/nobackup/nnsanche/NO_BHs/pioneer50h243.1536gst1bwK1/pioneer50h243.1536gst1bwK1.003456')
    s_cr = pynbody.load('/nobackup/ibutsky/simulations/patient0_new/pioneer50h243.1536gst1bwK1BH.002432')
    s_agncr = pynbody.load('/nobackup/ibutsky/simulations/patient0_agncr/pioneer50h243.1536gst1bwK1BH.%06d'%(output))
    plt.xlim(0, 12)
#    sfh_p0noBH = pp.sfh(s_sn, label = 'Thermal SNe only')
    sfh_cr = pp.sfh(s_cr, label = 'CR SNe feedback')
    sfh_p0 = pp.sfh(s_agn, label = 'P0')
    sfh_agncr = pp.sfh(s_agncr, label = 'P0+CR')
    
    plt.legend()
    plt.savefig('sfh_compare_%06d.png'%(output))
    print(sfh_cr[0])
#    print(sfh_cr.d)
    f = h5.File('sfh_data.h5', 'w')
    f.create_dataset("sfh_p0", data=sfh_p0[0])
    f.create_dataset("sfh_agncr", data=sfh_agncr[0])
    f.create_dataset("sfh_cr", data = sfh_cr[0])

    f.create_dataset("time_p0", data=sfh_p0[1])
    f.create_dataset("time_agncr", data=sfh_agncr[1])
    f.create_dataset("time_cr", data = sfh_cr[1])
    f.close()

def make_rotation_curves(s, output, cr):
    pynbody.analysis.angmom.faceon(s.g)
    if cr:
        s.g['eps'] = s.g['soft']
        s.d['eps'] = s.d['soft']
        s.s['eps'] = s.s['soft']

    p = pynbody.analysis.profile.Profile(s,min=.01,max=500,type='log',ndim=3)
    pg = pynbody.analysis.profile.Profile(s.g,min=.01,max=500,type='log',ndim=3)
    ps = pynbody.analysis.profile.Profile(s.s,min=.01,max=500,type='log',ndim=3)
    pd = pynbody.analysis.profile.Profile(s.d,min=.01,max=500,type='log',ndim=3)

    # make the plot
    plt.plot(p['rbins'],p['v_circ'],label='total')
    plt.plot(pg['rbins'],pg['v_circ'],label='gas')
    plt.plot(ps['rbins'],ps['v_circ'],label='stars')
    plt.plot(pd['rbins'],pd['v_circ'],label='dm')

    plt.xlabel('R [kpc]')
    plt.ylabel(r'$v_c$ [km/s]')
    plt.legend()
    if cr:
        plt.savefig('pioneer.%06d_rotation_curve.png'%(output))
    else:
        plt.savefig('pioneer.%06d_rotation_curve_nocr.png'%(output))
    plt.clf()

make_plots()
