import scipy.io
import numpy as np
import torch
import os

"""
POTENTIAL VALUES TO FILTER OUT:
Cc is sometimes 0.0000
IPAR is sometimes 0.0000

TODO: interpolate gap, gap-fills
"""

"""
The data is a dict from param_names to large arrays of values.
We transform it into a large array of dicts (and perform some translation and constant handling)
"""
def load_pipeline_data_dict(basepath, site_name, predictor_keys=None, constant_keys=None, output_keys=None, verbose=False): # NOTE: no real need to take 2 files into input, as the 2 files could be reconstructed from the site name
    
    # prepare paths
    predictor_path = os.path.join(basepath, f"Results_{site_name}.mat")
    observation_path = os.path.join(basepath, f"Res_{site_name}.mat")

    trunc = 100
    print_k_points=20
    
    # load dat from both mat files
    predictor_data = scipy.io.loadmat(predictor_path)
    observation_data = scipy.io.loadmat(observation_path)

    # build a dictionary with the desired values and the correct names
    ctx = "sun_H" # NOTE: Change here for different contexts (sunny/shaded, low/high vegetation)
    # sun_condition = if (LAI_L(i) > 0) && (Csno == 0) && (Cice == 0) TODO

    if predictor_keys is None:
        predictor_keys = ["Cc", "IPAR", "Csl", "ra", "rb", "Ts", "Pre", "Ds", "Psi_L", "Rn", "QG", "rs"] # NOTE: added rs for debugging 
    if constant_keys is None:
        constant_keys = ["Psi_sto_50", "Psi_sto_00", "CT", "Vmax", 
     "Ha", "FI", "Oa", "Do", "a1", "go", "gmes", "rjv", "DS"] # NOTE: missing DS. also why do we have integers for some constants?
    if output_keys is None:
        output_keys = ["LE_CORR"]


    # TODO: Very weird, understand why data points come by 2 all the time
    mapped_predictor_keys = key_mapping(predictor_keys, ctx)
    predictor_arrays = {
        k : torch.tensor(
            (predictor_data[k_m][:, 0] if predictor_data[k_m].shape[1]==2 else predictor_data[k_m])
            .astype(np.float32))
            .flatten()[:trunc] 
            for k,k_m in zip(predictor_keys, mapped_predictor_keys)}

    mapped_constant_keys = key_mapping(constant_keys, ctx)
    constants = {k : torch.tensor(predictor_data[k_m][0].astype(np.float32)).flatten()[:trunc] 
                for k,k_m in zip(constant_keys, mapped_constant_keys)}
    output_arrays = {k : torch.tensor(observation_data[k].astype(np.float32)).flatten()[:trunc] 
                    for k in output_keys}
    
    def valid(i):
        return (# True or ( # NOTE: this enables/disables the filtering
            not torch.isnan(output_arrays[output_keys[0]][i])  # Filter out nan output values! 
            and not predictor_arrays["Cc"][i] == 0
            and not predictor_arrays["IPAR"][i] == 0
        )

    # merge the dictionaries into an array of (input, output) tuples
    n = predictor_arrays[predictor_keys[0]].shape[0]
    n = trunc
    data = [
        ({ k : predictor_arrays[k][i] for k in predictor_arrays } | constants | {"DS_" : 0.5},
        output_arrays[output_keys[0]][i]) # For now, we assume a single output
        for i in range(n) if valid(i) # Filter out nan output values! 
    ]

    if verbose:
        #print(f"predictor_data={predictor_data}")
        #print(f"observation_data={observation_data}")
        #print(f"predictor_data[Ts]={predictor_data['Ta']}")
        #print(f"predictor_data[\"Ts\"][:{trunc}]{predictor_data['Ta'][:trunc]}")

        print(f" --------BASE SHAPES-------- ")
        for key in mapped_predictor_keys:
            print(f"Initial Shape for predictor {key}: {predictor_data[key].shape}")
        for key in output_keys:
            print(f"Initial Shape for output {key}: {observation_data[key].shape}")

        print(f"\n --------ARRAYS SHAPES-------- ")

        for k,v in predictor_arrays.items():
            print(f"predictor array [{k}] shape = {v.shape}")
            #print(f"predictor array [{k}] = {v}")
        for k,v in output_arrays.items():
            print(f"output array [{k}] shape = {v.shape}")
        #print(f"output_arrays={output_arrays}")

        print(f"n_points before filter: {len(predictor_arrays['Cc'])}")
        print(f"{print_k_points} first predictors Cc BEFORE: {predictor_arrays['Cc'][:print_k_points]}")
        print(f"{print_k_points} first output LE BEFORE: {output_arrays[output_keys[0]][:print_k_points]}")
        #print(f"first {print_k_points} predictors before filter: \n{predictor_arrays[:print_k_points]}")
        #print(f"first {print_k_points} outputs before filter: \n{output_arrays[:print_k_points]}")

        print("--------------FINAL FILTERED SHAPES-------------")

        print(f"n_points after filter: {len(data)}")
        print(f"first {print_k_points} data points:\n{data[:print_k_points]}")



    return data



def key_mapping(keys, ctx="sun_H"):
    if ctx == "sun_H":
        return [sun_h_map[k] for k in keys]
    elif ctx == "sun_L":
        return [sun_l_map[k] for k in keys]


sun_h_map = {
    "Cc":"Ci_sunH",
    "IPAR":"PAR_sun_H_final",
    "Csl":"Ca",
    "ra":"ran_H_final",
    "rb":"rb_H",
    "Ts":"Ta",
    "Pre":"Pre",
    "Ds":"Ds",
    "Psi_L":"Psi_l_H",
    "Psi_sto_50":"Psi_sto_50_H",
    "Psi_sto_00":"Psi_sto_00_H",
    "CT":"CT_H",
    "Vmax":"Vmax_H",
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
    "Psi_L":"Psi_l_L",
    "Psi_sto_50":"Psi_sto_50_L",
    "Psi_sto_00":"Psi_sto_00_L",
    "CT":"CT_L",
    "Vmax":"Vmax_L",
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
