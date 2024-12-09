import numpy as np
import torch
import os

from photosynthesis_biochemical import photosynthesis_biochemical
import train
import data

import test_differentiable_relations
import test_photosynthesis_biochemical

import plot


def debug_pipeline(base_path, site_name):
    output_keys=["LE"]
    LE_data = data.load_pipeline_data_dict(base_path, site_name, output_keys=output_keys) # NOTE: Could use LE_CORR
    pipeline = train.make_pipeline(None, None, output_rs=True)
    groundtruth_rs = []
    groundtruth_LE = []
    pred_rs = []
    pred_LE = []
    for i, (x, y) in enumerate(LE_data[:1000]):
        #print(f"---------------DATA POINT {i}----------------")
        x['CT'] = 3

        # test fix (to remove)
        #x['IPAR'] = 3.6639
        x['Vmax'] = 40.7293
        #x['Psi_L'] = -0.0229

        Q_LE, rs = pipeline(x)
        pred_rs.append(rs)
        pred_LE.append(Q_LE)
        groundtruth_rs.append(x['rs'])
        groundtruth_LE.append(y)
#
        #print(f"x = {x}")
        #print(f"pred rs={rs}")
        #print(f"pred Q_LE={Q_LE}")
        #print(f"obs rs = {x['rs']}")
        #print(f"obs Q_LE = {y}")

        #test_differentiable_relations.test_Q_LE_invertible([x])
        #test_differentiable_relations.test_Q_LE_using_inverted_PM([(x, y)])
        test_photosynthesis_biochemical.compare_rs([x])
    max_pred_Q_LE = np.array(pred_LE).max()
    max_groundtruth_Q_LE = np.array(groundtruth_LE).max()
    print(f"Max Pred Q_LE : {max_pred_Q_LE}")
    print(f"Max Groundtruth Q_LE : {max_groundtruth_Q_LE}")

def debug_data_loading(base_path, site_name):
    data.load_pipeline_data_dict(base_path, site_name, verbose=True) # NOTE: Could use LE_CORR

def plot_timeseries(base_path, site_name):
    ##### plot stuff
    print("Loading LE")
    LE_data = data.load_pipeline_data_dict(base_path, site_name, output_keys=["LE"])
    print("Loading LE_CORR")
    LE_corr_data = data.load_pipeline_data_dict(base_path, site_name, output_keys=["LE_CORR"])
    pipeline = train.make_pipeline(None, None)
    output_LE = []
    groundtruth_LE = []
    groundtruth_LE_CORR = []
    for x, y in LE_data:
        #print(f"x = {x}")
        x['CT'] = 3
        output_LE.append(pipeline(x))
        groundtruth_LE.append(y)
    for x, y in LE_corr_data:
        groundtruth_LE_CORR.append(y)
    plot.time_series_plot(timeseries=[(output_LE, "pred_LE"), (groundtruth_LE, "obs_LE"), (groundtruth_LE_CORR, "obs_LE_CORR")], name="LE_comparaison", show=True)
    #####

def train_pipeline(base_path, site_name):
    pipeline_data = data.load_pipeline_data_dict(base_path, site_name, output_keys=["LE_CORR"])
    gsCO2_model, Vmax_model = train.train_pipeline(pipeline_data)


def main():
    base_path = os.path.expanduser(
        "~/epfl/semester_project/databases/T_C_Input_and_Output_Pure_Physics/"
    )
    site_name = "CH-Dav"
    #site_name = "DE-Tha"

    debug_pipeline(base_path, site_name)
    #debug_data_loading(base_path, site_name)
    #plot_timeseries(base_path, site_name)






if __name__ == "__main__":
    main()
