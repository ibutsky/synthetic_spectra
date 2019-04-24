import numpy as np
import h5py as h5
import sys

sys.path.append('../plotting')
import plotting_tools as pt

def calculate_mean_err(data, err, use_log = False):

    # get rid of isnan error values
    mask = (~np.isnan(err)) 
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

    col = -9999.
    colerr = 0.
    vel = -9999.
    velerr = 0.
    bval = -9999.
    bvalerr = 0.
    flag = -9999.

    # first see if there are any good detections       
    detected = (flag_list == 1) & (veeper_col_list > -9999)
    sat      = (flag_list == 9) 
    uplim    = (flag_list == 5)
    if len(veeper_col_list[detected]) > 0:
        json_mask = detected & (json_col_list > -9999)
        aodm_mask = detected & (aodm_col_list > -9999)
        temp_col = np.append(veeper_col_list[detected], 
                                   [aodm_col_list[aodm_mask], 
                                   json_col_list[json_mask]])
        temp_col_err = np.append(veeper_colerr_list[detected], 
                                      [aodm_colerr_list[aodm_mask],
                                      json_colerr_list[json_mask]])
        col, colerr = calculate_mean_err(temp_col, temp_col_err, use_log = True)
            
            #inds    = veeper_colerr_list[detected].argsort()
            #col     = veeper_col_list[detected][inds][0]
            #colerr  = veeper_colerr_list[detected][inds][0]

#            vel     =      vel_list[detected][0]
#            velerr  =  vel_err_list[detected][0]
#            bval    =     bval_list[detected][0]
#            bvalerr = bval_err_list[detected][0]
        vel, velerr = calculate_mean_err(vel_list[detected], vel_err_list[detected])
        bval, bvalerr = calculate_mean_err(bval_list[detected], bval_err_list[detected])
        flag = 1

    elif len(veeper_col_list[sat])  > 0:
        inds = aodm_col_list[sat].argsort()
        col    =    aodm_col_list[sat][inds][0]
        col    =    np.min(np.append(aodm_col_list[sat], json_col_list[sat]))
        colerr =    0. #aodm_colerr_list[sat][inds][0]

        vel     =      vel_list[sat][0]
        velerr  =  vel_err_list[sat][0]
        bval    =     bval_list[sat][0]
        bvalerr = bval_err_list[sat][0]

        vel, velerr= calculate_mean_err(vel_list[sat], vel_err_list[sat])
        bval, bvalerr = calculate_mean_err(bval_list[sat], bval_err_list[sat])

        flag = 9

    elif len(veeper_col_list[uplim]) > 0:
        col    =    np.max(np.append(aodm_col_list[uplim], json_col_list[uplim]))
        colerr =    0; #aodm_colerr_list[uplim][inds][-1]

        vel     = -9999.
        velerr  = 0.
        bval    = -9999.
        bvalerr = 0.

        flag = 5
        

    return col, colerr, vel, velerr, bval, bvalerr, flag



all_data = h5.File('../../data/analyzed_spectra/combined_spectra.h5', 'r')
all_orientations, all_models, all_times, all_impacts, all_ions, veeper_cols, veeper_colerrs, \
    vels, velerrs, bvals, bvalerrs, aodm_cols, aodm_colerrs, json_cols, json_colerrs, aodm_flags = \
    pt.load_data(['orientation', 'model', 'time', 'impact', 'ion', 'col', 'colerr', 'vel', 'velerr', \
                  'bval', 'bvalerr', 'col_aodm', 'col_aodm_err', 'col_json', 'col_json_err', 'flag_aodm'], \
                     use_filtered = False)



orientation_list = [];      model_list = [];   time_list = [];  redshift_list = [];
impact_list      = [];        ion_list = [];    col_list = [];      flag_list = [];
sigcol_list      = [];       bval_list = [];    vel_list = [];   sigbval_list = [];
sigvel_list      = [];   ion_name_list = [];   


working_dir = '../../data/analyzed_spectra'
views = ['face', 'edge_theta0', 'edge_theta1.0', 'edge_theta1.5']
models = ['anisd', 'stream']
impacts = np.arange(10, 110, 10)
times = [11.2]
redshifts = [0.2]
ions = ['HI', 'OVI', 'CII', 'CIII', 'SiII', 'SiIII', 'SiIV',  'NIII', 'NV']
#ions = ['SiIII']
#models = ['stream']
#views = ['edge_theta1.0']
for view in views:
    for model in models:
        for impact in impacts:
            for time, redshift in zip(times, redshifts):
                for ion in ions:
                    mask = (all_orientations == view) & (all_models == model) & (all_impacts == impact) &\
                        (all_times == time) & (all_ions == ion) 
                    if len(all_ions[mask]) > 0:
                        col, colerr, vel, velerr, bval, bvalerr, use_flag = \
                            best_measurement(veeper_cols[mask], veeper_colerrs[mask], \
                            aodm_cols[mask], aodm_colerrs[mask], json_cols[mask], \
                            json_colerrs[mask], aodm_flags[mask], vels[mask], \
                            velerrs[mask], bvals[mask], bvalerrs[mask])
                      
                        if col > 0:
                            orientation_list.append( view)
                            model_list.append(      model)
                            time_list.append(        time)
                            redshift_list.append(redshift)
                            impact_list.append(    impact)
                            ion_name_list.append(     ion)

                            col_list.append(       col)
                            sigcol_list.append( colerr)
                            vel_list.append(       vel)
                            sigvel_list.append( velerr)
                            bval_list.append(     bval)
                            sigbval_list.append(bvalerr)
                            flag_list.append( use_flag)
                     
                    



spec_outfile = h5.File('filtered_spectra.h5', 'w')

dataset_names = ['impact', 'time', 'redshift', 'col', 'colerr', 'bval', 'bvalerr', 'vel', 'velerr', 'flag']
datasets = [impact_list, time_list, redshift_list, col_list, sigcol_list, bval_list, sigbval_list, vel_list, sigvel_list,flag_list]
# first save the numerical data   
for dset, data in zip(dataset_names, datasets):
    spec_outfile.create_dataset(dset, data = data)


# then save string-type data in a special way   
dt = h5.special_dtype(vlen=str)
dataset_names = ['model', 'orientation', 'ion']
datasets = [model_list, orientation_list, ion_name_list]
for dset, data in zip(dataset_names, datasets):
    current_dset = spec_outfile.create_dataset(dset, (len(data),), dtype=dt)
    for i in range(len(data)):
        current_dset[i] = data[i].replace(" ", "")

spec_outfile.close()
