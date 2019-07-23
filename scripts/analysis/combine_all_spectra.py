#### Reads in all of the ".dat" files generated by Veeper
#### and combines them into a single file, "all_spectra.h5"

import numpy as np
import glob
import os
import sys
import h5py as h5 

import eqwrange as eqw
import spectrum_analysis_tools as spa

master_ion_list = ['HI', 'OVI', 'CII', 'CIII', 'SiII', 'SiIII', 'SiIV', 'NV']

work_dir = '../../data/analyzed_spectra'
spec_outfile = h5.File('%s/combined_spectra.h5'%(work_dir), 'w')

model_list   = [];   redshift_list = [];   impact_list = [];     label_list = [];
col_list     = [];     sigcol_list = [];     bval_list = [];   sigbval_list = [];
vel_list     = [];     sigvel_list = [];      ion_list = [];    ray_id_list = [];
eqw_list     = [];     sigeqw_list = [];    lncol_list = [];  siglncol_list = []; 
velcent_list = [];   velwidth_list = [];   ewjson_list = []; sigewjson_list = []; 
coljson_list = []; sigcoljson_list = []; restwave_list = []; flag_aodm_list = [];


# go to the directory where the analyzed spectra reside
os.chdir(work_dir)
dummy = -9999.


spec_files = glob.glob('COS-FUV*')
for spec in spec_files:
    if not os.path.isdir(spec) or not spa.spec_ready_for_analysis(spec):
        print('Skipping %s\n'%(spec))
        continue

    model, redshift, impact, ray_id = spa.extract_spec_info(spec)

    veeper_fn = '%s/compiledVPoutputs.dat'%(spec)
    json_fn = '%s/%s_lineids.json'%(spec, spec)
    json_out = '%s/json_eqw.dat'%(spec)
    aodm_fn = '%s/%s_ibnorm.fits'%(spec, spec)
    aodm_plot_dir = '%s/aodm_plots'%(spec)
    if not os.path.isdir(aodm_plot_dir):
        os.mkdir(aodm_plot_dir)

    veeper_ions, veeper_restwaves, veeper_cols, veeper_colerr, veeper_bvals, \
        veeper_bvalerr, veeper_vels, veeper_velerr, veeper_label =  eqw.load_veeper_fit(veeper_fn)
    json_ions, json_restwaves, json_eqw, json_eqwerr, json_col, json_colerr = eqw.json_eqw(json_fn, aodm_fn, json_out)
    
    for i in range(len(master_ion_list)):
        ion = master_ion_list[i].replace(" ", "")
        all_restwaves = spa.all_restwaves(ion)
        for rw in all_restwaves:
            # assume the number of components is 1, to start
            num_comps = 1
            index = (veeper_ions == ion)# & (veeper_restwaves == rw)
            if ion in veeper_ions:# and veeper_ions[index].size > 0:
                if len(veeper_ions[index]) > 1:
                    num_comps += veeper_ions[index].size - 1
                restwave_list = np.append(restwave_list,   num_comps*[rw])
                ion_list      = np.append(ion_list,       num_comps*[ion])
                col_list      = np.append(col_list,         veeper_cols[index])
                sigcol_list   = np.append(sigcol_list,   veeper_colerr[index])
                bval_list     = np.append(bval_list,      veeper_bvals[index])
                sigbval_list  = np.append(sigbval_list, veeper_bvalerr[index])
                vel_list      = np.append(vel_list,         veeper_vels[index])
                sigvel_list   = np.append(sigvel_list,   veeper_velerr[index])
                label_list    = np.append(label_list,     veeper_label[index])

            else:
                restwave_list = np.append(restwave_list,    rw)
                ion_list      = np.append(ion_list,        ion)
                col_list      = np.append(col_list,      dummy)
                sigcol_list   = np.append(sigcol_list,   dummy)
                bval_list     = np.append(bval_list,     dummy)
                sigbval_list  = np.append(sigbval_list,  dummy)
                vel_list      = np.append(vel_list,      dummy)
                sigvel_list   = np.append(sigvel_list,   dummy)
                label_list    = np.append(label_list,     "--")

            json_index = (json_ions == ion) & (json_restwaves == rw)

            if ion in json_ions and len(json_ions[json_index]) > 0:
                if len(json_ions[json_index]) > num_comps:
                    ewjson_list     = np.append(ewjson_list,        json_eqw[json_index][:num_comps])
                    sigewjson_list  = np.append(sigewjson_list,  json_eqwerr[json_index][:num_comps])
                    coljson_list    = np.append(coljson_list,       json_col[json_index][:num_comps])
                    sigcoljson_list = np.append(sigcoljson_list, json_colerr[json_index][:num_comps])
                elif len(json_ions[json_index]) == 1 and num_comps > 1:
                    ewjson_list     = np.append(ewjson_list,        [json_eqw[json_index]]*num_comps)
                    sigewjson_list  = np.append(sigewjson_list,  [json_eqwerr[json_index]]*num_comps)
                    coljson_list    = np.append(coljson_list,       [json_col[json_index]]*num_comps)
                    sigcoljson_list = np.append(sigcoljson_list, [json_colerr[json_index]]*num_comps)
                else:
                    ewjson_list     = np.append(ewjson_list,        json_eqw[json_index])
                    sigewjson_list  = np.append(sigewjson_list,  json_eqwerr[json_index])
                    coljson_list    = np.append(coljson_list,       json_col[json_index])
                    sigcoljson_list = np.append(sigcoljson_list, json_colerr[json_index])
            else:
                ewjson_list     = np.append(ewjson_list,     num_comps*[dummy])
                sigewjson_list  = np.append(sigewjson_list,  num_comps*[dummy])
                coljson_list    = np.append(coljson_list,    num_comps*[dummy])
                sigcoljson_list = np.append(sigcoljson_list, num_comps*[dummy])

                
            impact_list      = np.append(impact_list,     num_comps*[impact])
            model_list       = np.append(model_list,      num_comps*[model])
            ray_id_list      = np.append(ray_id_list,     num_comps*[ray_id])
            redshift_list    = np.append(redshift_list,   num_comps*[redshift])

            eqws, sigeqw, lncol, siglncol, flag_aodm, velcent, velwidth = \
                eqw.find_ion_limits(ion, aodm_fn, restwave = rw, redshift = redshift,\
                                        silent = 1, plots = 0, plot_dir = aodm_plot_dir, \
                                        vrange = (-200, 200), sat_limit = 0.1)
            flag_aodm_list= np.append(flag_aodm_list, num_comps*[flag_aodm[0]])
            eqw_list      = np.append(eqw_list,            num_comps*[eqws[0]])
            sigeqw_list   = np.append(sigeqw_list,       num_comps*[sigeqw[0]])
            lncol_list    = np.append(lncol_list,         num_comps*[lncol[0]])
            siglncol_list = np.append(siglncol_list,   num_comps*[siglncol[0]])
            
print(label_list)
print(len(label_list), len(model_list))
dataset_names = ['impact', 'ray_id', 'redshift', 'restwave', 'col_veeper', 'col_err_veeper', 'bval', \
                  'bval_err', 'vel', 'vel_err', 'flag', 'eqw_aodm', 'eqw_err_aodm', 'col_aodm', \
                   'col_err_aodm', 'col_json', 'col_err_json', 'eqw_json', 'eqw_err_json']
datasets = [impact_list, ray_id_list, redshift_list, restwave_list, col_list, sigcol_list, bval_list, sigbval_list, vel_list, sigvel_list, \
                flag_aodm_list, eqw_list, sigeqw_list, lncol_list, siglncol_list, coljson_list, \
                sigcoljson_list, ewjson_list, sigewjson_list]

# first save the numerical data   
for dset, data in zip(dataset_names, datasets):
    print(dset)
    spec_outfile.create_dataset(dset, data = data)


# then save string-type data in a special way     
dt = h5.special_dtype(vlen=str)
dataset_names = ['model', 'ion', 'label']
datasets = [model_list, ion_list, label_list]
for dset, data in zip(dataset_names, datasets):
    print(dset)
    current_dset = spec_outfile.create_dataset(dset, (len(data),), dtype=dt)
    for i in range(len(data)):
        current_dset[i] = data[i].replace(" ", "")


spec_outfile.close()
            
                    
                    



