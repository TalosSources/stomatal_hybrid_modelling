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
    #output_keys=["LE_CORR"]
    output_keys=["LE"]
    LE_data = data.load_pipeline_data_dict(base_path, site_name, output_keys=output_keys) # NOTE: Could use LE_CORR
    pipeline = train.make_pipeline(None, None, output_rs=True)
    groundtruth_rs = []
    groundtruth_LE = []
    pred_rs = []
    pred_LE = []
    for i, data_point in enumerate(LE_data[:100]):
        print(f"---------------DATA POINT {i}----------------")

        # fix some ill values TODO: Do that elsewhere
        for ctx in data_point[0].keys():
            #print(f"fixing ctx {ctx}")
            data_point[0][ctx]['CT'] = 3

            # test fix (to remove)
            #x['IPAR'] = 3.6639
            #print(f"point before: {data_point[0][ctx]}")
            #data_point[0][ctx]['Vmax'] = torch.tensor([40.7293])
            #print(f"point after: {data_point[0][ctx]}")
            #data_point[0][ctx]['Psi_L'] = torch.tensor([-0.0229])


        #print(f"the datapoint is: {data_point}")

        #print(f"\n\nresulting values: \n")

        Q_LE = pipeline(data_point)
        # pred_rs.append(rs)
        pred_LE.append(Q_LE)
        # groundtruth_rs.append(x['rs'])
        groundtruth_LE.append(data_point[2])
#
        #print(f"x = {x}")
        #print(f"pred rs={rs}")
        #print(f"pred Q_LE={Q_LE}")
        #print(f"obs rs = {x['rs']}")
        #print(f"obs Q_LE = {data_point[2]}")

        #test_differentiable_relations.test_Q_LE_invertible([x])
        #test_differentiable_relations.test_Q_LE_using_inverted_PM([(x, y)])
        #test_photosynthesis_biochemical.compare_rs([x])
    max_pred_Q_LE = np.array(pred_LE).max()
    max_groundtruth_Q_LE = np.array(groundtruth_LE).max()
    print(f"Max Pred Q_LE : {max_pred_Q_LE}")
    print(f"Max Groundtruth Q_LE : {max_groundtruth_Q_LE}")
    min_pred_Q_LE = np.array(pred_LE).min()
    min_groundtruth_Q_LE = np.array(groundtruth_LE).min()
    print(f"Min Pred Q_LE : {min_pred_Q_LE}")
    print(f"Min Groundtruth Q_LE : {min_groundtruth_Q_LE}")
    plot.time_series_plot(timeseries=[(pred_LE, "pred_LE"), (groundtruth_LE, "obs_LE")], name="LE_output_input_comparaison", show=True)

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
    pipeline_data = data.load_pipeline_data_dict(base_path, site_name, output_keys=["LE"])
    gsCO2_model, Vmax_model = train.train_pipeline(pipeline_data)


def main():
    base_path = os.path.expanduser(
        "~/epfl/semester_project/databases/T_C_PIPELINE_DATA/"
    )
    site_name = "CH-Dav"
    #site_name = "DE-Tha"
    #site_name = "DE-Hai"
    #site_name = "US-Me2"

    #debug_pipeline(base_path, site_name)
    #debug_data_loading(base_path, site_name)
    #plot_timeseries(base_path, site_name)
    train_pipeline(base_path, site_name)





if __name__ == "__main__":
    main()
