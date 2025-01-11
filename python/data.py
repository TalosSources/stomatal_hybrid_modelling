import scipy.io
import numpy as np
import torch
import os

"""
This module offers methods to load T&C matlab files res.mat and results.mat,
and process them into handy and flexible dicts that can be used by the training pipeline.
"""

"""
POTENTIAL VALUES TO FILTER OUT:
Cc is sometimes 0.0000
IPAR is sometimes 0.0000 (and NaN)
Negative Q_LE is bad
in pb, An is sometimes 0, and sometimes even NaN [INVESTIGATED: CAUSED BY IPAR]
does every timestep have exactly 1 valid ctx? at least 1?

NOTE: can interpolate gap, gap-fills
"""


def load_pipeline_data_dict_from_all_sites(config, predictor_keys=None, constant_keys=None, 
                            output_keys=None, global_keys=None, verbose=False):

    # prepare paths
    base_path = os.path.expanduser(config.base_path)
    data = []
    for site_name in config.sites:
        predictor_path = os.path.join(base_path, f"Results_{site_name}.mat")
        observation_path = os.path.join(base_path, f"Res_{site_name}.mat")

        load_pipeline_data_dict_single_site(predictor_path, observation_path, data, config.nPoints, 
            predictor_keys=predictor_keys, constant_keys=constant_keys, output_keys=output_keys,
            global_keys=global_keys, verbose=verbose, exclude_negative_outputs=config.exclude_negative)
        
    return data


"""
The data is a dict from param_names to large arrays of values.
We transform it into a large array of dicts (and perform some translation and constant handling)
"""
def load_pipeline_data_dict_single_site(predictor_path, observation_path, data, nPoints = None, predictor_keys=None, 
                        constant_keys=None, output_keys=None, global_keys=None, verbose=False, exclude_negative_outputs=False):
    
    # load dat from both mat files
    predictor_data = scipy.io.loadmat(predictor_path)
    observation_data = scipy.io.loadmat(observation_path)

    # potentially truncate the site data
    site_data_length = observation_data["LE"].shape[0]
    if nPoints is None:
        nPoints = site_data_length
    else:
        nPoints = min(nPoints, site_data_length)

    # Prepare the default keys
    if predictor_keys is None:
        predictor_keys = ["Cc", "IPAR", "Csl", "ra", "rb", "Ts", "Pre", 
                          "Ds", "Psi_L", "Rn", "QG", "Vmax"] 
    if constant_keys is None:
        constant_keys = ["Psi_sto_50", "Psi_sto_00", "CT", 
     "Ha", "FI", "Oa", "Do", "a1", "go", "gmes", "rjv", "DS"]
    if output_keys is None:
        output_keys = ["LE"]
    if global_keys is None:
        global_keys = ["EIn_H", "EIn_L", "EG", "ELitter", "ESN", "ESN_In", 
                       "EWAT", "EICE", "EIn_urb", "EIn_rock"]  

    ctxs =  ["sun_H", "sun_L", "shd_H", "shd_L"]

    # Load the predictor arrays: For each ctx, for each pipeline predictor, a tensor
    predictor_arrays = {
        ctx: {
            k : torch.tensor(
                (predictor_data[k_m][:, 0:1] if predictor_data[k_m].shape[1]==2 else predictor_data[k_m])
                .astype(np.float32))
                .flatten()[:nPoints] 
                for k,k_m in zip(predictor_keys, key_mapping(predictor_keys, ctx))
        } for ctx in ctxs
    }

    # Load the constants array: For each ctx, for each constant, a single value
    constants = {
        ctx: {
            k : torch.tensor(predictor_data[k_m][0, 0:1].astype(np.float32))
                    for k,k_m in zip(constant_keys, key_mapping(constant_keys, ctx))
        } for ctx in ctxs
    }

    # Load the global variable arrays: For each global key, a tensor
    global_arrays = {
        k : torch.tensor(predictor_data[k].astype(np.float32)).flatten()[:nPoints]
        for k in global_keys
    }

    # Load the output arrays: For each output key, a tensor
    output_arrays = {
        k : torch.tensor(observation_data[k].astype(np.float32)).flatten()[:nPoints] 
        for k in output_keys
    }

    # Check whether a data point is valid.
    def is_valid_timestep(i):
        return  (
            not torch.isnan(output_arrays[output_keys[0]][i])  # Filter out nan Q_LE values! 
            and ((not exclude_negative_outputs) or (not (output_arrays[output_keys[0]][i] < 0.))) # Filter out negative Q_LE values! actually should we? our pipeline can actually output negative values
            # NOTE: Perhaps add some more checks if needed 
        )
    
    # Check whether a ctx for a specific timestep is valid.
    def is_valid_context(ctx, i):
        return (
            predictor_arrays[ctx]['ra'][i] != 0. # ra is in the denom of the PM equation: can't be 0.
            and not predictor_arrays[ctx]['IPAR'][i].isnan() # IPAR is used in the pb pipeline.
            # NOTE: Add other checks
            #and not predictor_arrays["Cc"][i] == 0 # Uncomment to filter out values where cc is 0 
            #and not predictor_arrays["IPAR"][i] == 0 # Uncomment to filter out values where IPAR is 0 
            and not constants[ctx]['CT'].isnan() # AU-TTE has nan values for L constants
        )
        
    

    # merge the dictionaries into an array of (input, output) tuples, convenient for the pipeline

    # Output and Global don't depend on the context
    # Predictors and constants do.
    # Data will have this shape:
    # An n-length array of pairs (x_tuple, y_output)
    # Each x_tuple is (single_pipeline_dict, global_dict)
    # global_dict maps global variables -> single number
    # single_pipeline_dict maps ctx -> predictor_dict
    # output is a single number (for now?)
    start_idx = len(data)
    for i in range(nPoints):
        # perform ctx checks, to find valid contexts
        valid_ctxs = [ctx for ctx in ctxs if is_valid_context(ctx, i)]

        # perform the global check
        is_valid = is_valid_timestep(i) and len(valid_ctxs) > 0

        # if the timestep is valid, include the valid contexts
        if is_valid:
            data.append((
                (
                    {ctx: 
                        { k : predictor_arrays[ctx][k][i].unsqueeze(0) for k in predictor_arrays[ctx] } | constants[ctx] | {"DS_" : 0.5}
                        for ctx in valid_ctxs
                    }, 
                    {k : global_arrays[k][i].unsqueeze(0) for k in global_arrays},
                ), 
                output_arrays[output_keys[0]][i], 
            ))



    if verbose:

        print_k_points=1

        site_name = predictor_path[-10:-4]
        print(f"---------LOADING DATA FOR SITE {site_name}---------")

        #print(f"predictor_data={predictor_data}")
        #print(f"observation_data={observation_data}")
        #print(f"predictor_data[Ts]={predictor_data['Ta']}")
        #print(f"predictor_data[\"Ts\"][:{trunc}]{predictor_data['Ta'][:trunc]}")
        #print(f"predictor_data[\"Ta\"] = {predictor_data['Ta'].shape}")
        #print(f"predictor_data[\"Psi_sto_50_H\"] = {predictor_data['Psi_sto_50_H'].shape}")

        #print(f" --------BASE SHAPES-------- ")
        #for key in ["Vmax_H", "Vmax_L"]:
        #    print(f"Initial Shape for predictor {key}: {predictor_data[key].shape}")
        #for key in output_keys:
        #    print(f"Initial Shape for output {key}: {observation_data[key].shape}")

        #print(f"\n --------ARRAYS SHAPES-------- ")

        # for ctx, preds in predictor_arrays.items():
        #     print(f"For context {ctx}:")
        #     for k,v in preds.items():
        #         print(f"    predictor array [{k}] shape = {v.shape}")
        #     #print(f"predictor array [{k}] = {v}")
        #for k,v in output_arrays.items():
        #    print(f"output array [{k}] shape = {v.shape}")
        #print(f"output_arrays={output_arrays}")

        #print(f"n_points before filter: {len(predictor_arrays['Cc'])}")
        #print(f"{print_k_points} first predictors Cc BEFORE: {predictor_arrays['Cc'][:print_k_points]}")
        #print(f"{print_k_points} first output LE BEFORE: {output_arrays[output_keys[0]][:print_k_points]}")
        #print(f"first {print_k_points} predictors before filter: \n{predictor_arrays[:print_k_points]}")
        #print(f"first {print_k_points} outputs before filter: \n{output_arrays[:print_k_points]}")

        #print("--------------FINAL FILTERED SHAPES-------------")

        print(f"found constants = {constants}")

        print(f"n_points after filter: {len(data) - start_idx}")
        print(f"first {print_k_points} data points:")
        for p in data[start_idx:start_idx+print_k_points]:
            print(f"{p}\n-----------------------------------\n")



def key_mapping(keys, ctx="sun_H"):
    if ctx == "sun_H":
        return [sun_h_map[k] for k in keys]
    elif ctx == "sun_L":
        return [sun_l_map[k] for k in keys]
    elif ctx == "shd_H":
        return [shd_h_map[k] for k in keys]
    elif ctx == "shd_L":
        return [shd_l_map[k] for k in keys]
    else:
        raise ValueError("key_mapping: Invalid context given")


sun_h_map = {
    "Cc":"Ci_sunH",
    "IPAR":"PAR_sun_H_final",
    "Csl":"Ca",
    "ra":"ran_H_final",
    "rb":"rb_H",
    "Ts":"Ta",
    "Pre":"Pre",
    "Ds":"Ds",
    "Psi_L":"Psi_sun_H_final",
    "Psi_sto_50":"Psi_sto_50_H",
    "Psi_sto_00":"Psi_sto_00_H",
    "CT":"CT_H",
    "Vmax":"Vmax_sun_H_final",
    "Ha":"Ha_H",
    "FI":"FI_H",
    "Oa":"Oa",
    "Do":"Do_H",
    "a1":"a1_H",
    "go":"go_H",
    "gmes":"gmes_H",
    "rjv":"rjv_H",
    "DS":"DSE_H", # NOTE: Hand-added
    "QG":"G", # NOTE: Hand-added
    "Rn":"Rn", # NOTE: Hand-added
    "rs" : "rs_sunH", # NOTE: Hand-added
}

sun_l_map = {
    "Cc":"Ci_sunL",
    "IPAR":"PAR_sun_L_final",
    "Csl":"Ca",
    "ra":"ran_L_final",
    "rb":"rb_L",
    "Ts":"Ta",
    "Pre":"Pre",
    "Ds":"Ds",
    "Psi_L":"Psi_sun_L_final",
    "Psi_sto_50":"Psi_sto_50_L",
    "Psi_sto_00":"Psi_sto_00_L",
    "CT":"CT_L",
    "Vmax":"Vmax_sun_L_final",
    "Ha":"Ha_L",
    "FI":"FI_L",
    "Oa":"Oa",
    "Do":"Do_L",
    "a1":"a1_L",
    "go":"go_L",
    "gmes":"gmes_L",
    "rjv":"rjv_L",
    "DS":"DSE_L", # NOTE: Hand-added
    "QG":"G", # NOTE: Hand-added
    "Rn":"Rn", # NOTE: Hand-added
    "rs" : "rs_sunL", # NOTE: Hand-added
}

shd_h_map = {
    "Cc":"Ci_shdH",
    "IPAR":"PAR_shd_H_final",
    "Csl":"Ca",
    "ra":"ran_H_final",
    "rb":"rb_H",
    "Ts":"Ta",
    "Pre":"Pre",
    "Ds":"Ds",
    "Psi_L":"Psi_sun_H_final",
    "Psi_sto_50":"Psi_sto_50_H",
    "Psi_sto_00":"Psi_sto_00_H",
    "CT":"CT_H",
    "Vmax":"Vmax_shd_H_final",
    "Ha":"Ha_H",
    "FI":"FI_H",
    "Oa":"Oa",
    "Do":"Do_H",
    "a1":"a1_H",
    "go":"go_H",
    "gmes":"gmes_H",
    "rjv":"rjv_H",
    "DS":"DSE_H", # NOTE: Hand-added
    "QG":"G", # NOTE: Hand-added
    "Rn":"Rn", # NOTE: Hand-added
    "rs" : "rs_shdH", # NOTE: Hand-added
}

shd_l_map = {
    "Cc":"Ci_shdL",
    "IPAR":"PAR_shd_L_final",
    "Csl":"Ca",
    "ra":"ran_L_final",
    "rb":"rb_L",
    "Ts":"Ta",
    "Pre":"Pre",
    "Ds":"Ds",
    "Psi_L":"Psi_sun_L_final",
    "Psi_sto_50":"Psi_sto_50_L",
    "Psi_sto_00":"Psi_sto_00_L",
    "CT":"CT_L",
    "Vmax":"Vmax_shd_L_final",
    "Ha":"Ha_L",
    "FI":"FI_L",
    "Oa":"Oa",
    "Do":"Do_L",
    "a1":"a1_L",
    "go":"go_L",
    "gmes":"gmes_L",
    "rjv":"rjv_L",
    "DS":"DSE_L", # NOTE: Hand-added
    "QG":"G", # NOTE: Hand-added
    "Rn":"Rn", # NOTE: Hand-added
    "rs" : "rs_shdL", # NOTE: Hand-added
}
