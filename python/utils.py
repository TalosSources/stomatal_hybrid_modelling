from itertools import product

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
