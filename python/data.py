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

TODO: interpolate gap, gap-fills
"""

"""
The data is a dict from param_names to large arrays of values.
We transform it into a large array of dicts (and perform some translation and constant handling)
"""
def load_pipeline_data_dict(config, predictor_keys=None, constant_keys=None, 
                            output_keys=None, global_keys=None, verbose=False):
    
    # prepare paths
    base_path = os.path.expanduser(config.base_path)
    site_name = config.sites[0] # NOTE: Could add the option to merge data from multiple sites
    predictor_path = os.path.join(base_path, f"Results_{site_name}.mat")
    observation_path = os.path.join(base_path, f"Res_{site_name}.mat")

    # load dat from both mat files
    predictor_data = scipy.io.loadmat(predictor_path)
    observation_data = scipy.io.loadmat(observation_path)

    # potentially truncate the site data
    nPoints = config.nPoints
    if nPoints is None:
        observation_data["LE"].shape[0]

    # Prepare the default keys
    if predictor_keys is None:
        predictor_keys = ["Cc", "IPAR", "Csl", "ra", "rb", "Ts", "Pre", 
                          "Ds", "Psi_L", "Rn", "QG", "rs", "Vmax"] # NOTE: added rs for debugging 
    if constant_keys is None:
        constant_keys = ["Psi_sto_50", "Psi_sto_00", "CT", 
     "Ha", "FI", "Oa", "Do", "a1", "go", "gmes", "rjv", "DS"]
    if output_keys is None:
        output_keys = ["LE_CORR"]
    if global_keys is None:
        global_keys = ["EIn_H", "EIn_L", "EG", "ELitter", "ESN", "ESN_In", 
                       "EWAT", "EICE", "EIn_urb", "EIn_rock"]  

    # sun_condition = if (LAI_L(i) > 0) && (Csno == 0) && (Cice == 0) TODO
    # TODO: Data points come by 2 because there's 2 types of high vegetation in some sites. 
    # In this case, either we ignore the site, or we can treat them as 2 separate contexts?
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

    # Check whether a data point is valid. maybe TO CONFIG which check to perform
    def is_valid_timestep(i):
        return  (
            not torch.isnan(output_arrays[output_keys[0]][i])  # Filter out nan Q_LE values! 
            #and not output_arrays[output_keys[0]][i] < 0. # Filter out negative Q_LE values! actually should we? our pipeline can actually output negative values TO CONFIG
            # NOTE: Perhaps add some more checks if needed 
        )
    
    # Check whether a ctx for a specific timestep is valid. maybe TO CONFIG which check to perform
    def is_valid_context(ctx, i):
        return (
            predictor_arrays[ctx]['ra'][i] != 0. # ra is in the denom of the PM equation: can't be 0.
            and not predictor_arrays[ctx]['IPAR'][i].isnan() # IPAR is used in the pb pipeline.
            # NOTE: Add other checks
            #and not predictor_arrays["Cc"][i] == 0 # Filter out values where cc is 0 (for which context???) I guess, it should 0 for all contexts?
            #and not predictor_arrays["IPAR"][i] == 0 # Filter out values where IPAR is 0 (for which context???)
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
    data = []
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

        #print(f"predictor_data={predictor_data}")
        #print(f"observation_data={observation_data}")
        #print(f"predictor_data[Ts]={predictor_data['Ta']}")
        #print(f"predictor_data[\"Ts\"][:{trunc}]{predictor_data['Ta'][:trunc]}")
        print(f"predictor_data[\"Ta\"] = {predictor_data['Ta'].shape}")
        print(f"predictor_data[\"Psi_sto_50_H\"] = {predictor_data['Psi_sto_50_H'].shape}")

        print(f" --------BASE SHAPES-------- ")
        for key in ["Vmax_H", "Vmax_L"]:
            print(f"Initial Shape for predictor {key}: {predictor_data[key].shape}")
        for key in output_keys:
            print(f"Initial Shape for output {key}: {observation_data[key].shape}")

        print(f"\n --------ARRAYS SHAPES-------- ")

        for ctx, preds in predictor_arrays.items():
            print(f"For context {ctx}:")
            for k,v in preds.items():
                print(f"    predictor array [{k}] shape = {v.shape}")
            #print(f"predictor array [{k}] = {v}")
        for k,v in output_arrays.items():
            print(f"output array [{k}] shape = {v.shape}")
        #print(f"output_arrays={output_arrays}")

        #print(f"n_points before filter: {len(predictor_arrays['Cc'])}")
        #print(f"{print_k_points} first predictors Cc BEFORE: {predictor_arrays['Cc'][:print_k_points]}")
        #print(f"{print_k_points} first output LE BEFORE: {output_arrays[output_keys[0]][:print_k_points]}")
        #print(f"first {print_k_points} predictors before filter: \n{predictor_arrays[:print_k_points]}")
        #print(f"first {print_k_points} outputs before filter: \n{output_arrays[:print_k_points]}")

        print("--------------FINAL FILTERED SHAPES-------------")

        print(f"n_points after filter: {len(data)}")
        print(f"first {print_k_points} data points:\n{data[:print_k_points]}")



    return data



def key_mapping(keys, ctx="sun_H"):
    # NOTE: We could factorize this ugly code
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
