import os
import glob
import numpy as np
import sys

def fix_joebvp_file(fn, redshift):
    restwaves, cols, bvals, vels, nflag, bflag, vflag, vlim1, vlim2, wobs1, wobs2, z_comp = \
    np.loadtxt(fn, unpack=True, skiprows = 1, \
                usecols = (1,3,4,5,6,7,8,9, 10, 11, 12, 13), delimiter = '|')
    file_name, veeper_ions, rely, comment = \
    np.loadtxt(fn, unpack=True, skiprows = 1, usecols = (0, 14, 15, 16), dtype = 'str', delimiter = '|')

    if not isinstance(restwaves, np.ndarray):
        restwaves  = [restwaves]; cols = [cols]; bvals = [bvals];   vels = [vels]; nflag = [nflag]; 
        bflag      = [bflag];   vflag = [vflag]; vlim1 = [vlim1]; vlim2 = [vlim2]; wobs1 = [wobs1]; wobs2 = [wobs2]; 
        z_comp     = [z_comp]; file_name = [file_name]; veeper_ions = [veeper_ions]; 
        rely = [rely]; comment = [comment];

    zsys = redshift*np.ones(len(restwaves))    
    output_fn = open(fn, 'w')
    output_fn.write('specfile|restwave|zsys|col|bval|vel|nflag|bflag|vflag|vlim1|vlim2|wobs1|wobs2|z_comp|trans|rely|comment\n')
    for i in range(len(restwaves)):
        output_fn.write('%s|%f|%f|%f|%f|%f|%i|%i|%i|%f|%f|%f|%f|%f|%s|%s|%s\n'%\
                    (file_name[i], restwaves[i], zsys[i], cols[i], bvals[i], vels[i], nflag[i], bflag[i], \
                     vflag[i], vlim1[i], vlim2[i], wobs1[i], wobs2[i], z_comp[i], veeper_ions[i], rely[i], comment[i]))
    output_fn.close()


os.chdir('../../data/analyzed_spectra/')

#views =  ['edge_theta1.5']
#models = ['anisd']
#impacts = [50]

views = ['face', 'edge_theta0', 'edge_theta1.0', 'edge_theta1.5']
models = ['anisd', 'stream']
impacts = np.arange(10, 210, 10)

views = ['edge_theta1.0']
models = ['anisd']
impacts = [40]
time = 11.2
z = 0.2

for view in views:
    for model in models:
        for impact in impacts:
            folder_name = 'COS-FUV_%s_%s_%0.1fGyr_r%ikpc/'%(view, model, time, impact)
            print(folder_name)
            if not os.path.isdir(folder_name):
                continue
            os.chdir(folder_name)
            if os.path.isfile('compiledVPoutputs.dat'):
                os.remove('compiledVPoutputs.dat')
            joebvp_files = glob.glob('*.joebvp')
            for fn in joebvp_files:
                print(fn)
                sys.stdout.flush()
                fix_joebvp_file(fn, z)
            os.chdir('../')
os.chdir('../../scripts/analysis/')
