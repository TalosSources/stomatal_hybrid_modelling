import numpy as np
import torch
import os

from photosynthesis_biochemical import photosynthesis_biochemical
import train
import data

import plot


def debug_pipeline(base_path, site_name):
    output_keys=["LE"]
    LE_data = data.load_pipeline_data_dict(base_path, site_name, output_keys=output_keys) # NOTE: Could use LE_CORR
    pipeline = train.make_pipeline(None, None, output_rs=True)
    groundtruth_rs = []
    groundtruth_LE = []
    pred_rs = []
    pred_LE = []
    for x, y in LE_data[:10]:
        x['CT'] = 3
        Q_LE, rs = pipeline(x)
        pred_rs.append(rs)
        pred_LE.append(Q_LE)
        groundtruth_rs.append(x['rs'])
        groundtruth_LE.append(y)

        print(f"x = {x}")
        print(f"pred rs={rs}")
        print(f"pred Q_LE={Q_LE}")
        print(f"obs rs = {x['rs']}")
        print(f"obs Q_LE = {y}")


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






if __name__ == "__main__":
    main()
