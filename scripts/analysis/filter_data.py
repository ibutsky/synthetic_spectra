import numpy as np
import h5py as h5
import sys

sys.path.append('../plotting')
import plotting_tools as pt

def calculate_mean_err(data, err, use_log = False):

    # get rid of isnan error values
#    mask = (~np.isnan(err)) 
    mask = err > 0
    err = err[mask]
    data = data[mask]

    if use_log:
        data = 10**data
        err = 10**err
    mean = np.mean(data)
    std  = np.sqrt( np.sum( err**2  ) / len(err) )
    if use_log:
        mean = np.log10(mean)
        std = np.log10(std)
    return mean, std

def best_measurement(veeper_col_list, veeper_colerr_list, aodm_col_list, aodm_colerr_list, \
                         json_col_list, json_colerr_list, flag_list, \
                         vel_list, vel_err_list, bval_list, bval_err_list):

    col = [-9999.]
    colerr = [0.]
    vel = [-9999.]
    velerr = [0.]
    bval = [-9999.]
    bvalerr = [0.]
    flag = [-9999.]

    # first see if there are any good detections       
    detected = (flag_list == 1) & (veeper_col_list > -9999)
    sat      = (flag_list == 9) 
    uplim    = (flag_list == 5) | (veeper_col_list == -9999)
    if len(veeper_col_list[detected]) > 0:
        col    = veeper_col_list[detected]
        colerr = veeper_colerr_list[detected]
        vel    = vel_list[detected]
        velerr = vel_err_list[detected]
        bval   = bval_list[detected]
        bvalerr = bval_err_list[detected]
        flag   = flag_list[detected]

#
#        json_mask = detected & (json_col_list > -9999)
#        aodm_mask = detected & (aodm_col_list > -9999)
#        temp_col = np.append(veeper_col_list[detected], 
#                                   [aodm_col_list[aodm_mask], 
#                                   json_col_list[json_mask]])
#        temp_col_err = np.append(veeper_colerr_list[detected], 
#                                      [aodm_colerr_list[aodm_mask],
#                                      json_colerr_list[json_mask]])
#        col, colerr = calculate_mean_err(temp_col, temp_col_err, use_log = True)
            
            #inds    = veeper_colerr_list[detected].argsort()
            #col     = veeper_col_list[detected][inds][0]
            #colerr  = veeper_colerr_list[detected][inds][0]

#            vel     =      vel_list[detected][0]
#            velerr  =  vel_err_list[detected][0]
#            bval    =     bval_list[detected][0]
#            bvalerr = bval_err_list[detected][0]
#        vel, velerr = calculate_mean_err(vel_list[detected], vel_err_list[detected])
#        bval, bvalerr = calculate_mean_err(bval_list[detected], bval_err_list[detected])
#        flag = 1

    elif len(veeper_col_list[sat])  > 0:
        inds = aodm_col_list[sat].argsort()
        col    =    aodm_col_list[sat][inds][0]
        col    =    [np.min(np.append(aodm_col_list[sat], json_col_list[sat]))]
        colerr =    [0.] #aodm_colerr_list[sat][inds][0]

        vel     =      [vel_list[sat][0]]
        velerr  =  [vel_err_list[sat][0]]
        bval    =    [ bval_list[sat][0]]
        bvalerr = [bval_err_list[sat][0]]

        #vel, velerr= calculate_mean_err(vel_list[sat], vel_err_list[sat])
        #bval, bvalerr = calculate_mean_err(bval_list[sat], bval_err_list[sat])

        flag = [9]

    elif len(veeper_col_list[uplim]) > 0:
        col    =    [np.max(np.append(aodm_col_list[uplim], json_col_list[uplim]))]
        colerr =    [0.] #aodm_colerr_list[uplim][inds][-1]

        vel     =[ -9999.]
        velerr  = [0.]
        bval    = [-9999.]
        bvalerr = [0.]

        flag = [5]
        

    return col, colerr, vel, velerr, bval, bvalerr, flag



all_data = h5.File('../../data/analyzed_spectra/combined_spectra.h5', 'r')
all_models, all_ray_ids, all_redshifts, all_impacts, all_ions, veeper_cols, veeper_colerrs, \
    vels, velerrs, bvals, bvalerrs, aodm_cols, aodm_colerrs, json_cols, json_colerrs, aodm_flags = \
    pt.load_data(['model', 'ray_id', 'redshift', 'impact', 'ion', 'col_veeper', 'col_err_veeper', 'vel', 'vel_err', \
                  'bval', 'bval_err', 'col_aodm', 'col_err_aodm', 'col_json', 'col_err_json', 'flag'], \
                     use_filtered = False)



ray_id_list = [];           model_list = [];   redshift_list = [];  
impact_list      = [];        ion_list = [];    col_list = [];      flag_list = [];
sigcol_list      = [];       bval_list = [];    vel_list = [];   sigbval_list = [];
sigvel_list      = [];   ion_name_list = [];   


working_dir = '../../data/analyzed_spectra'
models = ['P0', 'tempest']
ray_ids = np.arange(20)
redshifts = [0.25]
ions = ['HI', 'OVI', 'CII', 'CIII', 'CIV', 'SiII', 'SiIII', 'SiIV', 'NIII', 'NV']

for model in models:
    for redshift in redshifts:
        for ray_id in ray_ids:
            for ion in ions:
                mask = (all_models == model) & (all_ray_ids == ray_id) &\
                        (all_redshifts == redshift) & (all_ions == ion) 
                if len(all_ions[mask]) > 0:
                    col, colerr, vel, velerr, bval, bvalerr, use_flag = \
                            best_measurement(veeper_cols[mask], veeper_colerrs[mask], \
                            aodm_cols[mask], aodm_colerrs[mask], json_cols[mask], \
                            json_colerrs[mask], aodm_flags[mask], vels[mask], \
                            velerrs[mask], bvals[mask], bvalerrs[mask])
                      
                    if col[0] > 0:
                        col_list      = np.append(col_list, col)
                        sigcol_list   = np.append(sigcol_list, colerr)
                        vel_list      = np.append(vel_list, vel)
                        sigvel_list   = np.append(sigvel_list, velerr)
                        bval_list     = np.append(bval_list, bval)
                        sigbval_list  = np.append(sigbval_list, bvalerr)
                        flag_list     = np.append(flag_list, use_flag)
                     
                        num = len(col)
                    
                        model_list       = np.append(model_list,  num*[model])
                        ray_id_list      = np.append(ray_id_list, num*[ray_id])
                        redshift_list    = np.append(redshift_list,  num*[redshift])

                        impact = all_impacts[mask][0]
                        impact_list      = np.append(impact_list,  num*[impact])
                        ion_name_list    = np.append(ion_name_list,  num*[ion])                    


print(len(impact_list), len(col_list), len(vel_list), len(model_list))
spec_outfile = h5.File('../../data/analyzed_spectra/filtered_spectra.h5', 'w')

dataset_names = ['impact', 'redshift', 'col', 'col_err', 'bval', 'bval_err', 'vel', 'vel_err', 'flag']
datasets = [impact_list, redshift_list, col_list, sigcol_list, bval_list, sigbval_list, vel_list, sigvel_list,flag_list]
# first save the numerical data   
for dset, data in zip(dataset_names, datasets):
    data = data.astype('float64')
    spec_outfile.create_dataset(dset, data = data)


# then save string-type data in a special way   
dt = h5.special_dtype(vlen=str)
dataset_names = ['model', 'ion']
datasets = [model_list, ion_name_list]
for dset, data in zip(dataset_names, datasets):
    current_dset = spec_outfile.create_dataset(dset, (len(data),), dtype=dt)
    for i in range(len(data)):
        current_dset[i] = data[i].replace(" ", "")

spec_outfile.close()
