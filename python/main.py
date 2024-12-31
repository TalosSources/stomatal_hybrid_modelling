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
    pipeline_data = data.load_pipeline_data_dict_from_all_sites(config.data, output_keys=["LE"], verbose=True)

    if config.train.hp_tuning:
        best_hps = train.perform_hp_tuning(config, pipeline_data)
    else:
        rs_model, _ = train.train_and_evaluate_pipeline(config, pipeline_data)

        # Save the rs_model
        torch.save(rs_model.state_dict(), config.model_weights_path)



def main():

    #config = OmegaConf.load("configs/default.yaml") # NOTE: Could be passed as arg, or maybe we can have a list of experiments to perform
    config = OmegaConf.load("configs/hyperparameter_tuning.yaml")
    train_and_eval_pipeline(config)





if __name__ == "__main__":
    main()
