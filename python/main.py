import numpy as np
import torch
import os

from photosynthesis_biochemical import photosynthesis_biochemical
import train
import data

import plot

import eval


def train_and_eval_pipeline(base_path, site_name):
    nPoints = 1000 # TO CONFIG -> and decide if we want to truncate, what to do in training?
    pipeline_data = data.load_pipeline_data_dict(base_path, site_name, output_keys=["LE"], nPoints=nPoints, verbose=True)
    gsCO2_model, Vmax_model = train.train_pipeline(pipeline_data)
    # TODO: Store the models?


def main():
    
    base_path = os.path.expanduser( # TO CONFIG
        "~/epfl/semester_project/databases/T_C_PIPELINE_DATA/"
    )

    # Choose a(some) site(s) # TO CONFIG
    site_name = "CH-Dav"
    #site_name = "DE-Tha"
    #site_name = "DE-Hai"
    #site_name = "US-Me2"

    train_and_eval_pipeline(base_path, site_name)





if __name__ == "__main__":
    main()
