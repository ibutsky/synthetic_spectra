import os
import glob
import numpy as np
import sys

def fix_joebvp_file(fn, redshift):
    restwaves, cols, bvals, nflag, bflag, vflag, vlim1, vlim2, z_comp = \
    np.loadtxt(fn, unpack=True, skiprows = 1, \
                usecols = (1,3,4, 6,7,8,9, 10, 13), delimiter = '|')
    file_name, veeper_ions, rely, comment = \
    np.loadtxt(fn, unpack=True, skiprows = 1, usecols = (0, 14, 15, 16), dtype = 'str', delimiter = '|')

    if not isinstance(restwaves, np.ndarray):
        restwaves  = [restwaves]; cols = [cols]; bvals = [bvals]; nflag = [nflag]; 
        bflag      = [bflag];   vflag = [vflag]; vlim1 = [vlim1]; vlim2 = [vlim2]; 
        z_comp     = [z_comp]; file_name = [file_name]; veeper_ions = [veeper_ions]; 
        rely = [rely]; comment = [comment];

    zsys = redshift*np.ones(len(restwaves))    
    vels = (z_comp - zsys) * 299792.458
    wobs1 = restwaves * (1 + z_comp + vlim1 / 299792.458)
    wobs2 = restwaves * (1 + z_comp + vlim2 / 299792.458)
    
    output_fn = open(fn, 'w')
    output_fn.write('specfile|restwave|zsys|col|bval|vel|nflag|bflag|vflag|vlim1|vlim2|wobs1|wobs2|z_comp|trans|rely|comment\n')
    for i in range(len(restwaves)):
        output_fn.write('%s|%f|%f|%f|%f|%f|%i|%i|%i|%f|%f|%f|%f|%f|%s|%s|%s\n'%\
                    (file_name[i], restwaves[i], zsys[i], cols[i], bvals[i], vels[i], nflag[i], bflag[i], \
                     vflag[i], vlim1[i], vlim2[i], wobs1[i], wobs2[i], z_comp[i], veeper_ions[i], rely[i], comment[i]))
    output_fn.close()



redshift = float(sys.argv[1])
joebvp_files = glob.glob("*.joebvp")

for fn in joebvp_files:
    fix_joebvp_file(fn, redshift)
