import os
import glob
import numpy as np
import sys
from astropy import constants as const

# speed of light in km/s
c = const.c.to('km/s').value

def fix_joebvp_file(fn, redshift):
    # if there's an empty joebvp file, remove the file and exit
    if sum(1 for line in open(fn)) == 1:
       os.remove(fn)
       return

    restwaves, cols, bvals, nflag, bflag, vflag, vlim1, vlim2, wobs1, wobs2, z_comp = \
    np.loadtxt(fn, unpack=True, skiprows = 1, \
                usecols = (1,3,4, 6,7,8,9,10,11,12,13), delimiter = '|')

    file_name, veeper_ions, rely, comment = \
    np.loadtxt(fn, unpack=True, skiprows = 1, usecols = (0, 14, 15, 16), dtype = 'str', delimiter = '|')

    if not isinstance(restwaves, np.ndarray):
        restwaves  = np.array([restwaves]); cols = np.array([cols]); bvals = np.array([bvals]); 
        nflag = np.array([nflag]);  bflag = np.array([bflag]); vflag = np.array([vflag]); 
        vlim1 = np.array([vlim1]);  vlim2 = np.array([vlim2]); z_comp = np.array([z_comp]); 
        wobs1 = np.array([wobs1]); wobs2 = np.array([wobs2]);
        file_name = np.array([file_name]); veeper_ions = np.array([veeper_ions]); 
        rely = np.array([rely]); comment = np.array([comment]);

    zsys = redshift*np.ones(len(restwaves))    
    vels = c * ((1 + z_comp) / (1 + zsys) - 1)
    
    dv = (wobs2 - wobs1)/restwaves/(1 + z_comp) * c / 2.0
    vlim1 = vels - dv
    vlim2 = vels + dv
    
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
