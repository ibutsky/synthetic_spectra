import os
import glob
import numpy as np
import sys
import shutil

def clean_compiled_VPoutputs(fn):
    ''' Removes repeat data entries in compiledVPoutputs.dat for readability '''

    restwaves, zsys, cols, sigcols, bvals, sigbvals, vels, sigvels,  nflag, bflag, vflag, \
        vlim1, vlim2, wobs1, wobs2, pix1, pix2, z_comp  =  np.loadtxt(fn, unpack=True, skiprows = 1, \
                usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18), delimiter = '|')
    file_name, veeper_ions, rely, comment = \
    np.loadtxt(fn, unpack=True, skiprows = 1, usecols = (0, 19, 20, 21), dtype = 'str', delimiter = '|')

    mask = sigcols > 0
    restwaves = restwaves[mask]; zsys = zsys[mask]; cols = cols[mask]; sigcols=sigcols[mask];
    bvals = bvals[mask]; sigbvals = sigbvals[mask]; vels = vels[mask]; sigvels=sigvels[mask];
    nflag = nflag[mask]; bflag = bflag[mask]; vflag = vflag[mask]; vlim1 = vlim1[mask];
    vlim2 = vlim2[mask]; wobs1 = wobs1[mask]; wobs2 = wobs2[mask]; pix1 = pix1[mask]; 
    pix2 = pix2[mask]; z_comp = z_comp[mask]; file_name = file_name[mask]; 
    veeper_ions = veeper_ions[mask]; rely = rely[mask]; comment = comment[mask];


    output_fn = open(fn, 'w')
    output_fn.write('specfile|restwave|zsys|col|sigcol|bval|sigbval|vel|sigvel|nflag|bflag|vflag|vlim1|vlim2|wobs1|wobs2|pix1|pix2|z_comp|trans|rely|comment\n')
    for i in range(len(restwaves)):
        output_fn.write('%s|%8.3f|%0.2f|%6.3f|%6.3f|%6.3f|%6.3f|%6.3f|%6.3f|%i|%i|%i|%0.1f|%0.1f|%0.6f|%0.6f|%5i|%5i|%0.6f|%5s|%s|%s\n'%\
                    (file_name[i], restwaves[i], zsys[i], cols[i], sigcols[i], bvals[i], sigbvals[i], vels[i], sigvels[i], \
                     nflag[i], bflag[i], vflag[i], vlim1[i], vlim2[i], wobs1[i], wobs2[i], pix1[i], pix2[i], z_comp[i], 
                     veeper_ions[i], rely[i], comment[i]))
    output_fn.close()


spec_files = glob.glob('../../data/analyzed_spectra/COS-FUV*')
#for spec in spec_files:
#    fn = '%s/compiledVPoutputs.dat'%(spec)
#    original_fn = '%s/originalVPoutputs.dat'%(spec)
#    if not os.path.isfile(original_fn):
#        shutil.copy(fn, original_fn)
#    clean_compiled_VPoutputs(fn)

