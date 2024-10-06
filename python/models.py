"""
Pytorch network are stored in this file
"""

import torch
import numpy as np

class FCN(torch.nn.Module):
    def __init__(self, layers, activation = torch.tanh):
        super(FCN, self).__init__()
        
        self.layers = torch.nn.ModuleList(
            torch.nn.Linear(layers[i], layers[i+1])
            for i in range(len(layers)-1)
        )

        print(f"Initialized layers: {self.layers}")
        
        self.activation = activation
        
    
    def forward(self, x):
        print(f"calling forward, x={x}")
        for i, layer in enumerate(self.layers):
            x = layer(x)
            print(f"using layer {layer}, getting x={x}")
            if i < len(self.layers) - 1: # NOTE: Last layer is linear, necessary to express arbitrary real valued functions
                x = self.activation(x)
                print(f"using activation, getting x={x}")
            
        print(f"returning {x}")
        return x

def gsCO2_model():
    input_dim = 6
    output_dim = 1
    return FCN([input_dim, 64, 64, 64, output_dim], torch.nn.Tanh()) # NOTE: relatively simple network, subject to change (activation?)