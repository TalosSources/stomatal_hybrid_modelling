from itertools import product
import numpy as np
import torch
import os
from omegaconf import OmegaConf

import models


def printGradInfo(x, name):
    print(f"{name} value: {x}")
    print(f"{name} type: {type(x)}")
    print(f"{name} requireGrad: {x.requires_grad}")
    print(f"{name} grad: {x.grad}")


"""
Given a dict param_name (str) -> list of possible values,
returns an iterator of all possible combinations,
i.e. an iterator of dicts param_name (str) -> single value
"""


def grid_search_iterator(possible_values):
    keys = possible_values.keys()
    value_combinations = product(*possible_values.values())

    for combination in value_combinations:
        yield dict(zip(keys, combination))


"""
Give a unique id to each possible combination of ctxs
"""


def compute_ctx_set_id(ctx_set):
    base_ctxs = ["sun_H", "sun_L", "shd_H", "shd_L"]
    id = 0
    for i, ctx in enumerate(base_ctxs):
        if ctx in ctx_set:
            id += 2**i
    return id


def group_by_ctx_id(data):
    # first use a dict because we don't which ctx sets are present
    grouped_data = {}
    for t in data:
        x, _ = t
        ctx_preds, _ = x
        id = compute_ctx_set_id(ctx_preds.keys())
        for p in ctx_preds.values():
            if p["CT"] == 3:  # differentiate between CT==3 and CT==4
                id += 2**4
            break
        if id not in grouped_data:
            grouped_data[id] = []
        grouped_data[id].append(t)

    # then build an ordered list with the map
    grouped_data = list(grouped_data.values())

    probs = []
    total_length = len(data)
    for i in range(len(grouped_data)):
        size = len(grouped_data[i])
        grouped_data[i] = np.array(grouped_data[i], dtype=object)
        probs.append(size / total_length)

    return grouped_data, probs


"""
given a list of datapoints in the usual format (tuple preds, output, with global, ctx, etc.) 
Returns a batch in the usual format, containing tensors instead of single points.
Assumes all given datapoints share the same context set
"""


def make_batch(data):
    # predictors is an array of size batch_size of (ctx_dict, global_dict) values
    # y is an array of size batch_size of y values
    predictors, outputs = zip(*data)

    # outputs_batch is already a y tensor in the form we want it
    outputs_tensor = torch.tensor(outputs)

    # ctx_dicts is an array of size batch_size of ctx_dict values, same for global_dicts
    ctx_dicts, global_dicts = zip(*predictors)

    # for each key present in global_dicts, we collect all the batch_size values for this key in a single tensor
    # and map the key to that tensor
    global_data = {
        key: torch.tensor([global_dict[key] for global_dict in global_dicts])
        for key in global_dicts[0].keys()
    }

    # for each ctx present in ctx_dicts, and for each key present for this ctx,
    # aggregate all these values into a tensor, and map the corresponding key to this tensor,
    # and the corresponding ctx to this predictor_dict.
    ctx_data = {
        ctx: {
            key: torch.tensor([ctx_dict[ctx][key] for ctx_dict in ctx_dicts])
            for key in ctx_dicts[0][ctx].keys()
        }
        for ctx in ctx_dicts[0].keys()
    }

    return ((ctx_data, global_data), outputs_tensor)


def make_eval_batches(data):
    grouped_data, _ = group_by_ctx_id(data)
    return [make_batch(ctx_data) for ctx_data in grouped_data]


def hp_set_to_string(hp_set):
    # hp set is param_name -> value.
    # the value is usually a number but not always?
    s = ""
    for param, val in hp_set.items():
        s += f"{param}-{val}_"
    return s


def save_losses(losses, results_path, suffix=""):
    losses = np.array(losses)
    if suffix != "":
        suffix = "_" + suffix
    np.save(os.path.join(results_path, f"losses{suffix}"), losses)


def save_multiple_benchmarks(benchmarks, results_path):
    path = os.path.join(results_path, "benchmark.txt")
    with open(path, "w") as text_file:
        for name, val in benchmarks.items():
            text_file.write(f"{name} : {val} [s]\n")


def save_single_benchmark(benchmark, results_path):
    path = os.path.join(results_path, "benchmark.txt")
    with open(path, "w") as text_file:
        text_file.write(f"{benchmark} [s]\n")


def save_coeffs_of_determination(coeffs, results_path):
    path = os.path.join(results_path, "coefficients_of_determination.txt")
    with open(path, "w") as text_file:
        for name, val in coeffs.items():
            text_file.write(f"{name} : {val}\n")


def save_best_hps(hps, results_path):
    path = os.path.join(results_path, "best_hyperparameters.txt")
    with open(path, "w") as text_file:
        for name, val in hps.items():
            text_file.write(f"{name} : {val}\n")


def load_model(config_name):
    # load the model weights
    results_path = os.path.join("results", config_name)
    rs_model_state_dict = torch.load(os.path.join(results_path, "model_weights.pt"))
    cfg = OmegaConf.load(f"configs/{config_name}.yaml")
    rs_model = models.rs_model_from_config(cfg.model)
    rs_model.load_state_dict(rs_model_state_dict)
    return rs_model


def trace_model(model, model_name):
    # put model in inference mode
    model.eval()

    # move to CPU
    model.to("cpu")

    # forward the model with dummy input
    X = torch.rand(1, 7)
    traced_model = torch.jit.trace(model.forward, X)
    traced_model.save(f"traced_models/traced_{model_name}.pt")
