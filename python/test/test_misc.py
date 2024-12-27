import train
import data
import plot

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