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

def gsCO2_model():
    input_dim = 6
    output_dim = 1
    return FCN([input_dim, 32, 32, output_dim], torch.nn.Tanh()) # NOTE: relatively simple network, subject to change (activation?)