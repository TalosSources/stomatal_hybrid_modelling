import scipy.io
import numpy as np
import torch
import os

def load_pipeline_data(basepath, site_name): # NOTE: no real need to take 2 files into input, as the 2 files could be reconstructed from the site name
    
    # prepare paths
    predictor_path = os.path.join(basepath, f"Results_{site_name}.mat")
    observation_path = os.path.join(basepath, f"Res_{site_name}.mat")
    
    # load dat from both mat files
    predictor_data = scipy.io.loadmat(predictor_path)
    observation_data = scipy.io.loadmat(observation_path)

    # build a dictionary with the desired values and the correct names
    ctx = "sun_H" # NOTE: Change here for different contexts (sunny/shaded, low/high vegetation)
    predictor_keys = ["Cc", "IPAR", "Csl", "ra", "rb", "Ts", "Pre", "Ds", "Psi_L"]
    constant_keys = ["Psi_sto_50", "Psi_sto_00", "CT", "Vmax", 
     "Ha", "FI", "Oa", "Do", "a1", "go", "gmes", "rjv"] # NOTE: missing DS. also why do we have integers for some constants?
    output_keys = ["LE"]

    mapped_predictor_keys = key_mapping(predictor_keys, ctx)
    predictor_arrays = {k : torch.tensor(predictor_data[k_m].astype(np.float32)).flatten() for k,k_m in zip(predictor_keys, mapped_predictor_keys)}
    mapped_constant_keys = key_mapping(constant_keys, ctx)
    constants = {k : torch.tensor(predictor_data[k_m].astype(np.float32)).flatten() for k,k_m in zip(constant_keys, mapped_constant_keys)}
    output_arrays = {k : torch.tensor(observation_data[k].astype(np.float32)) for k in output_keys}
    print(f"output_arrays={output_arrays}")

    # merge the dictionaries into an array of (input, output) tuples
    n = predictor_arrays[predictor_keys[0]].shape[0]
    data = [
        ({ k : predictor_arrays[k][i] for k in predictor_arrays } | constants |{"DS" : 0.5},
        output_arrays[output_keys[0]][i]) # For now, we assume a single output
        for i in range(n)
    ]

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
}
