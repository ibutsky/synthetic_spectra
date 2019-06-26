import os
import numpy as np
from joebvp import VPmeasure

def veeper_fit(view, model, time, r, working_dir = '../../data/analyzed_spectra/'):
    current_dir = os.getcwd()
    os.chdir(working_dir)
    basename = 'COS-FUV_%s_%s_%0.1fGyr_r%ikpc'%(view, model, time, r)
    if os.path.isdir(basename):
        os.chdir(basename)
        if os.path.isfile('compiledVPoutputs.dat'):
            print('Veeper batch fit already run for %s'%(basename))
        else:
            hack = False
            num_lines = sum(1 for line in open('flist'))
            if num_lines > 0:
                if num_lines == 1:
                    hack = True
                VPmeasure.batch_fit('%s_ibnorm.fits'%(basename), \
                                    'flist', hack=hack)
        os.chdir('../')
    else:
        print('No folder: %s'%(basename))
    os.chdir(current_dir)


view =  'edge_theta1.0'
#view = 'face'
model = 'anisd'
time = 11.2
impacts = np.arange(10,  110, 10)
impacts = [40]
#for impact in impacts:
 #   veeper_fit(view, model, time, impact)

views = ['face', 'edge_theta0', 'edge_theta1.0', 'edge_theta1.5']
models = ['anisd', 'stream']
impacts = np.arange(10, 210, 10)
time = 11.2

views = ['edge_theta1.0']
models = ['anisd']
impacts = [40]

temp = open('temp.dat', 'w')
for view in views:
    for model in models:
        for impact in impacts:
            temp.write('%s, %s, %i\n'%(view, model, impact))
            temp.flush()
            veeper_fit(view, model, time, impact)
