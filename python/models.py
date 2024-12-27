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
            
        print(f"Initialized layers: {self.layers}")
        
    
    def forward(self, x):
        to_print = f"calling forward, x={x}\n"
        for i, layer in enumerate(self.layers):
            x = layer(x)
            to_print += f"using layer {i}={layer}, getting x={x}\n"
            #to_print += f"using layer {layer} (matrix={layer.weight}), getting x={x}\n" # with layer debug
        to_print += f"returning {x}\n"
        #print(to_print)
        return x


def gsCO2_model():
    """
    Ideas for improvement:
    * weight decay? normalisation? that kind of stuff
    * random forest?
    * proper recurrent network?
    * read the ML course again
    """
    input_dim = 6
    output_dim = 1

    # Different possible architectures TO CONFIG probably, if I want to give this option
    #return FCN([input_dim, 128, 64, 32, output_dim], torch.nn.ReLU()) # NOTE: relatively simple network, subject to change (activation?)
    #return FCN([input_dim, 32, output_dim], torch.nn.ReLU()) # NOTE: relatively simple network, subject to change (activation?)
    #return FCN([input_dim, 32, 32, output_dim], torch.nn.ReLU()) # NOTE: relatively simple network, subject to change (activation?)
    return FCN([input_dim, 32, 32, output_dim], torch.nn.ReLU(), batch_norm=True) # NOTE: relatively simple network, subject to change (activation?)

    # Minimal test: Linear model
    #return FCN([input_dim, output_dim], torch.nn.ReLU())

def vm_model():
    input_dim = 2
    output_dim = 1
    return FCN([input_dim, 64, 64, 64, output_dim], torch.nn.ReLU())

