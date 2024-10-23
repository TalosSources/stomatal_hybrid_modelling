import scipy.io
import os

def load_pipeline_data(basepath, site_name): # NOTE: no real need to take 2 files into input, as the 2 files could be reconstructed from the site name
    
    predictor_path = os.path.join(basepath, f"Results_{site_name}.mat")
    observation_path = os.path.join(basepath, f"Res_{site_name}.mat")
    
    predictor_data = scipy.io.loadmat(predictor_path)
    observation_data = scipy.io.loadmat(observation_path)
    
    print(f"PredictorData = {predictor_data}")
    
    return data