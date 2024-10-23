import numpy as np
import torch
import os

from photosynthesis_biochemical import photosynthesis_biochemical
import train
import data





#model = train.train_gsco2()
#model = train.train_vm()

def main():
    base_path = "~/epfl/semester_project/databases/T_C_Input_and_Output_Pure_Physics/"
    base_path = os.path.expanduser(base_path)
    site_name = "CH-Dav"

    pipeline_data = data.load_pipeline_data(base_path, site_name)
    #gsCO2_model, Vmax_model = train.train_pipeline(data)


if __name__ == "__main__":
    main()
