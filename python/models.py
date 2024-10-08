"""
Pytorch network are stored in this file
"""

import torch
import numpy as np

class FCN(torch.nn.Module):
    def __init__(self, layers, activation = torch.tanh):
        super(FCN, self).__init__()
        
        self.layers = torch.nn.ModuleList(
            torch.nn.Linear(layers[i], layers[i+1], bias=True)
            for i in range(len(layers)-1)
        )

        print(f"Initialized layers: {self.layers}")
        
        self.activation = activation
        
    
    def forward(self, x):
        to_print = f"calling forward, x={x}\n"
        for i, layer in enumerate(self.layers):
            x = layer(x)
            to_print += f"using layer {layer}, getting x={x}\n"
            #to_print += f"using layer {layer} (matrix={layer.weight}), getting x={x}\n" # with layer debug
            if i < len(self.layers) - 1: # NOTE: Last layer is linear, necessary to express arbitrary real valued functions
                x = self.activation(x)
                to_print += f"using activation, getting x={x}\n"
        
        to_print += f"returning {x}\n"
        print(to_print)
        return x
    
"""
Class for simple creation of simple fully connected recurrent networks
"""
class FCRN(torch.nn.Module):
    def __init__(self, input_dim, output_dim, n_hidden=2, dim_hidden=32, activation = torch.tanh):
        super(FCRN, self).__init__()

        input_layer = [torch.nn.Linear(input_dim, dim_hidden, bias=True)]
        hidden_layers = [torch.nn.Linear(dim_hidden, dim_hidden, bias=True) for i in range(n_hidden)]
        output_layer = [torch.nn.Linear(dim_hidden, output_dim, bias=True)]
        
        self.layers = torch.nn.ModuleList(input_layer + hidden_layers + output_layer)

        print(f"Initialized layers: {self.layers}")
        
        self.activation = activation
        
    
    def forward(self, x):
        to_print = f"calling forward, x={x}\n"

        # INPUT LAYER
        x = self.activation(self.layers[0](x)) # NOTE: Could remove this line by starting with x = 0, and do everything in the loop
        
        # RECURRENT LAYERS
        for i in range(len(self.layers)-2):
            residue = self.activation(self.layers[i+1](x))
            x = x + residue # RECURRENT STEP
            to_print += f"using layer {self.layers[i+1]}, getting x={x}\n"
            #to_print += f"using layer {layer} (matrix={layer.weight}), getting x={x}\n" # with layer debug
        
        # OUTPUT LAYER
        x = self.layers[-1](x)

        to_print += f"returning {x}\n"
        print(to_print)
        return x

def gsCO2_model():
    input_dim = 6
    output_dim = 1
    return FCN([input_dim, 64, 64, 64, 64, output_dim], torch.nn.ReLU()) # NOTE: relatively simple network, subject to change (activation?)
    #return FCRN(input_dim, output_dim, 3, 64, torch.nn.Tanh())