import numpy as np
import torch
import os
import sys
from omegaconf import OmegaConf

from photosynthesis_biochemical import photosynthesis_biochemical
import train
import data

import plot
import utils

import eval


def train_and_eval_pipeline(config):
    pipeline_data = data.load_pipeline_data_dict_from_all_sites(config.data, output_keys=["LE"], verbose=True)

    results_path = os.path.join("results", config.name)
    os.makedirs(results_path, exist_ok=True)

    if config.train.hp_tuning:
        best_hps, losses, benchmarks, coeffs_determination = train.perform_hp_tuning(config, pipeline_data)
        for hp_set_string, loss in losses.items():
            utils.save_losses(loss, results_path, suffix=hp_set_string)
        utils.save_multiple_benchmarks(benchmarks, results_path)
        utils.save_coeffs_of_determination(coeffs_determination, results_path)
        utils.save_best_hps(best_hps, results_path)
    else:
        rs_model, losses, benchmark = train.train_and_evaluate_pipeline(config, pipeline_data)
        utils.save_losses(losses, results_path)
        utils.save_single_benchmark(benchmark, results_path)

        # Save the rs_model
        torch.save(rs_model.state_dict(), os.path.join(results_path, "model_weights.pt"))



def main():

    if len(sys.argv) < 2:
        config_path = 'configs/best_model.yaml'
    else:
        config_path = sys.argv[1]

    config = OmegaConf.load(config_path)
    train_and_eval_pipeline(config)



if __name__ == "__main__":
    main()
