from itertools import product
import numpy as np

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
    base_ctxs = ['sun_H', 'sun_L', 'shd_H', 'shd_L']
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
        if id not in grouped_data:
            grouped_data[id] = []
        grouped_data[id].append(t)

    # then build an ordered list with the map
    grouped_data = list(grouped_data.values()) 

    probs = []
    total_length = len(data)
    for i in range(len(grouped_data)):
        size = len(grouped_data[i])
        grouped_data[i] = np.array(grouped_data[i])
        probs.append(size / total_length)

    return grouped_data, probs