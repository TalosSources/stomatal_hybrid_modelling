"""
This file contains method useful for evaluating, interpreting and plotting 
the results and components of the experiments.
This is used in writing the report for this project.
"""
import numpy as np
import os
from omegaconf import OmegaConf

import torchviz

import plot
import pipelines
import models
import data

def generate_hp_tuning_loss_plot():
    results_path = "results/lr_wd_tuning"
    filenames = [
        #"losses_lr-0.01_weight_decay-0_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        #"losses_lr-0.001_weight_decay-0_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        #"losses_lr-0.0001_weight_decay-0_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        "losses_lr-0.01_weight_decay-0.0001_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        "losses_lr-0.001_weight_decay-0.0001_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        "losses_lr-0.0001_weight_decay-0.0001_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        #"losses_lr-0.01_weight_decay-1e-06_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        #"losses_lr-0.001_weight_decay-1e-06_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        #"losses_lr-0.0001_weight_decay-1e-06_n_hidden-4_hidden_size-128_batch_norm-True_.npy"
    ]

    losses = [
        np.load(os.path.join(results_path, filename)) for filename in filenames
    ]

    labels = [
        #"lr=1e-2, lambda=0",
        #"lr=1e-3, lambda=0",
        #"lr=1e-4, lambda=0",
        "lr=1e-2, lambda=1e-4",
        "lr=1e-3, lambda=1e-4",
        "lr=1e-4, lambda=1e-4",
        #"lr=1e-2, lambda=1e-6",
        #"lr=1e-3, lambda=1e-6",
        #"lr=1e-4, lambda=1e-6",
    ]

    coefficients = [
        0.25,
        0.36,
        0.78
    ]

    plot.plot_losses(losses, labels, coefficients=coefficients, minmax=False, smoothing=0.01)

def generate_torch_viz_graph():
    # TODO: For now it only shows the model architecture, and some long line.
    # First, decide if I want to keep the model architecture. Perhaps not, it's not the focus
    # Then, give some of the predictors to torchviz, to display the actual Q_LE PM pipeline
    # If possible, find a datapoint with 2 or 4 valid contexts. (further? perhaps the crazy AU-Stp?)
    # Then, try to remove that long line of ADDs and stuff.
    # Finally, try to give better names to the nodes if possible

    cfg = OmegaConf.load("configs/default.yaml")
    cfg.data.nPoints = 10
    cfg.data.sites = ["CH-Dav"]
    cfg.model.batch_norm = False
    rs_model = models.rs_model_from_config(cfg.model)
    pipeline = pipelines.make_pipeline(rs_model, None, None)
    datapoints = data.load_pipeline_data_dict_from_all_sites(cfg.data)
    datapoint = datapoints[0]
    x, y = datapoint
    output = pipeline(x)
    graph = torchviz.make_dot(output, params=dict(rs_model.named_parameters()))
    graph.render("torchviz", view=True)

#generate_hp_tuning_loss_plot()
generate_torch_viz_graph()