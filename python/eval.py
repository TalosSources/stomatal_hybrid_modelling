import torch
import numpy as np
from sklearn.metrics import r2_score
from tqdm import tqdm

import utils

"""
General model evaluation function.
args:
    eval_data -> iterable of (x, y) tuples
    model -> function mapping x to a predicted y_hat
    criterion -> loss function taking y and y_hat as input
"""
def eval_general(model, eval_data, error_function):
    total_error = 0
    n_points = len(eval_data)
    eval_data = utils.make_eval_batches(eval_data)
    with torch.no_grad():
        for x, y in eval_data:
            n_points_batch = len(y)
            output = model(x)
            error = error_function(y, output)
            total_error += (n_points_batch/n_points) * error
    
    return total_error

def coefficient_of_determination(model, test_data):
    test_data = utils.make_eval_batches(test_data)
    ys = torch.tensor([])
    y_preds = torch.tensor([])
    for x, y in test_data:
        y_pred = model(x)
        # ys.append(y)
        ys = torch.cat([ys, y])
        # y_preds.append(y_pred)
        y_preds = torch.cat([y_preds, y_pred])

    return r2_score(ys.detach().numpy(), y_preds.detach().numpy())