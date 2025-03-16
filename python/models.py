"""
Pytorch network architectures are stored in this file
"""

import torch

class FCN(torch.nn.Module):
    def __init__(self, layers, activation = torch.tanh, batch_norm=False):
        super(FCN, self).__init__()

        module_list = []
        for i in range(len(layers) - 1):
            module_list.append(torch.nn.Linear(layers[i], layers[i+1], bias=True))
            if i < len(layers) - 2: # don't apply activation and batchNorm after the last linear layer
                if batch_norm:
                    module_list.append(torch.nn.BatchNorm1d(layers[i+1]))
                module_list.append(activation)
        self.layers = torch.nn.ModuleList(module_list)
            
    def forward(self, x):
        for i, layer in enumerate(self.layers):
            x = layer(x)
        return x


def rs_model_from_config(config):
    return rs_model(config.hidden_size, config.n_hidden, config.activation, config.batch_norm)

def rs_model(hidden_size, n_hidden, activation, batch_norm):
    input_dim = 7
    output_dim = 1

    # Architecture is specified in the config
    activation_map = { # Can be extended
        "ReLU" : torch.nn.ReLU
    }
    layers = [input_dim] + [hidden_size]*n_hidden + [output_dim]
    return FCN(layers, activation_map[activation](), batch_norm=batch_norm)

def gsCO2_model_from_config(config):
    return gsCO2_model(config.hidden_size, config.n_hidden, config.activation, config.batch_norm)

def gsCO2_model(hidden_size, n_hidden, activation, batch_norm):
    """
    Ideas for improvement:
    * weight decay? normalisation? that kind of stuff
    * random forest?
    * proper recurrent network?
    * read the ML course again
    """
    input_dim = 6 # the only difference with the rs_model
    output_dim = 1

    # Architecture is specified in the config
    activation_map = {
        "ReLU" : torch.nn.ReLU
    }
    layers = [input_dim] + [hidden_size]*n_hidden + [output_dim]
    return FCN(layers, activation_map[activation](), batch_norm=batch_norm)

