import os
import numpy as np
import glob
from joebvp import VPmeasure

def single_veeper_fit(basename, joebvp_list = 'flist',  working_dir = '../../data/analyzed_spectra', force_override = False):
    current_dir = os.getcwd()
    os.chdir('%s/%s'%(working_dir, basename))
    if os.path.isfile('compiledVPoutputs.dat') and not force_override:
        print('Veeper batch fit already run for %s'%(basename))
    else:
        VPmeasure.batch_fit('%s_ibnorm.fits'%(basename),joebvp_list)
    os.chdir(current_dir)


def all_veeper_fit(working_dir = '../../data/analyzed_spectra', joebvp_list = 'flist', force_override = False):
    current_dir = os.getcwd()
    os.chdir(working_dir)
    all_files = glob.glob('*')
    for basename in all_files:
        print(basename)
        if os.path.isdir(basename) and os.path.isfile('%s/%s'%(basename,joebvp_list)):
            if os.stat('%s/%s'%(basename, joebvp_list)).st_size == 0:
                print("WARNING: %s is empty; Skipping veeper fit for %s\n"%(joebvp_list, basename))
            else:
                single_veeper_fit(basename, joebvp_list = joebvp_list, force_override = force_override)
    os.chdir(current_dir)

#basename = 'COS-FUV_P0_z0.17_1'
#single_veeper_fit(basename, force_override = True)

all_veeper_fit()
