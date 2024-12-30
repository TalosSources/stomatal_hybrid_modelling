import numpy as np
import torch
import os
from omegaconf import OmegaConf

from photosynthesis_biochemical import photosynthesis_biochemical
import train
import data

import plot

import eval


def train_and_eval_pipeline(config):
    pipeline_data = data.load_pipeline_data_dict(config.data, output_keys=["LE"], verbose=True)

    if config.train.hp_tuning:
        best_hps = train.perform_hp_tuning(config, pipeline_data)
    else:
        gsCO2_model, Vmax_model = train.train_and_evaluate_pipeline(config, pipeline_data)

    # TODO: Store the models?



def main():

    config = OmegaConf.load("configs/default.yaml") # NOTE: Could be passed as arg, or maybe we can have a list of experiments to perform
    train_and_eval_pipeline(config)





if __name__ == "__main__":
    main()
