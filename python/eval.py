import torch
import numpy as np
from sklearn.metrics import r2_score

"""
General model evaluation function.
args:
    eval_data -> iterable of (x, y) tuples
    model -> function mapping x to a predicted y_hat
    criterion -> loss function taking y and y_hat as input
"""
def eval_general(model, eval_data, error_function):
    total_error = 0
    with torch.no_grad():
        for x, y in eval_data:
            output = model(x)
            #print(f"evaluating with x={x}, obtained out={output}")
            error = error_function(y, output)
            total_error += error
    
    total_error /= len(eval_data) # NOTE: Optional, computes the average loss

    return total_error

def coefficient_of_determination(model, test_data):
    ys = []
    y_preds = []
    for x, y in test_data:
        y_pred = model(x)
        ys.append(y)
        y_preds.append(y_pred)

    return r2_score(np.array(ys), np.array(y_preds))