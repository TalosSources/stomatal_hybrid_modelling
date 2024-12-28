import torch

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