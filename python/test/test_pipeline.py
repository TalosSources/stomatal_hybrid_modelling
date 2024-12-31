import numpy as np
import os

from data import load_pipeline_data_dict_single_site
from pipelines import make_pipeline
from plot import time_series_plot

# Take inspiration of those
def debug_pipeline(base_path, site_name):
    #output_keys=["LE_CORR"]
    output_keys=["LE"]
    LE_data = load_pipeline_data_dict_single_site(base_path, site_name, output_keys=output_keys, nPoints=10000) # NOTE: Could use LE_CORR
    pipeline = make_pipeline(None, None, output_rs=False)
    groundtruth_rs = []
    groundtruth_LE = []
    pred_rs = []
    pred_LE = []
    for i, data_point in enumerate(LE_data[:10000]):
        print(f"---------------DATA POINT {i}----------------")

        #print(f"the datapoint is: {data_point}")

        #print(f"\n\nresulting values: \n")
        x, y = data_point

        Q_LE = pipeline(x)
        # pred_rs.append(rs)
        pred_LE.append(Q_LE)
        # groundtruth_rs.append(x['rs'])
        groundtruth_LE.append(y)

        #print(f"x = {x}")
        #print(f"pred rs={rs}")
        #print(f"pred Q_LE={Q_LE}")
        #print(f"obs rs = {x['rs']}")
        #print(f"obs Q_LE = {data_point[2]}")

        #test_differentiable_relations.test_Q_LE_invertible([x])
        #test_differentiable_relations.test_Q_LE_using_inverted_PM([(x, y)])
        #test_photosynthesis_biochemical.compare_rs([x])
    #max_pred_Q_LE = np.array(pred_LE).max()
    #max_groundtruth_Q_LE = np.array(groundtruth_LE).max()
    #print(f"Max Pred Q_LE : {max_pred_Q_LE}")
    #print(f"Max Groundtruth Q_LE : {max_groundtruth_Q_LE}")
    #min_pred_Q_LE = np.array(pred_LE).min()
    #min_groundtruth_Q_LE = np.array(groundtruth_LE).min()
    #print(f"Min Pred Q_LE : {min_pred_Q_LE}")
    #print(f"Min Groundtruth Q_LE : {min_groundtruth_Q_LE}")
    print(f"------------------------------------")
    time_series_plot(timeseries=[(pred_LE, "pred_LE"), (groundtruth_LE, "obs_LE")], name="LE_output_input_comparaison", show=True)


base_path = os.path.expanduser("~/epfl/semester_project/databases/T_C_PIPELINE_DATA/")
site_name = "CH-Dav"
debug_pipeline(base_path, site_name)