import numpy as np
import torch
import os

from photosynthesis_biochemical import photosynthesis_biochemical
import train
import data


def main():
    base_path = os.path.expanduser(
        "~/epfl/semester_project/databases/T_C_Input_and_Output_Pure_Physics/"
    )
    site_name = "CH-Dav"

    pipeline_data = data.load_pipeline_data_dict(base_path, site_name)
#    pipeline_data = data.load_pipeline_data_tensor(base_path, site_name)
    gsCO2_model, Vmax_model = train.train_pipeline(pipeline_data)


if __name__ == "__main__":
    main()
