import torch

"""
General model evaluation function.
args:
    eval_data -> iterable of (x, y) tuples
    model -> function mapping x to a predicted y_hat
    criterion -> loss function taking y and y_hat as input
"""
def eval_general(model, eval_data, criterion):
    total_loss = 0
    with torch.no_grad():
        for x, y in eval_data:
            output, _ = model(x)
            loss = criterion(y, output)
            total_loss += loss
    
    total_loss /= len(eval_data) # NOTE: Optional, computes the average loss

    return total_loss