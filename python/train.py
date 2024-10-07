from photosynthesis_biochemical import photosynthesis_biochemical

import models
import torch
import numpy as np

def train():

    # training parameters
    epochs = ...

    # tensor containing the tunable parameters of the semi-empirical model
    parameters = ...

    # the loss used to compare ground-truth and predictions
    criterion = ...

    # gradient descent optimizer that adjusts the parameters
    optimizer = ...

    for i in range(epochs):

        # load the arrays corresponding to function parameters and ground-truth output values
        # (for example, some of the parameters of photosynthesis_biochemical, and gsCO2 values)
        x = ...
        y = ...


        Cc,IPAR,Csl,ra,rb,Ts,Pre,Ds, Psi_L,Psi_sto_50,Psi_sto_00, CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,gmes,rjv = x

        result = photosynthesis_biochemical(Cc,IPAR,Csl,ra,rb,Ts,Pre,Ds, Psi_L,Psi_sto_50,Psi_sto_00, CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,gmes,rjv)
        CcF,An,rs,Rdark,F755nm,GAM,gsCO2 = result

        loss = criterion(gsCO2, y)
        loss.backward()
        optimizer.step()
        optimizer.zerograd()

        ... # and so on, use some sensible gradient descent optimizer with some sensible loss


def train_general(
    model,
    criterion,
    optimizer,
    epochs,
    data_iterator,
    batch_size=1
):

    losses : list(float) = []

    for i in range(epochs):
        
        optimizer.zero_grad()
        
        x_batch = []
        y_batch = []
        for _ in range(batch_size):
            x, y = next(data_iterator)
            x_batch.append(x)
            y_batch.append(y)
        x = torch.tensor(x_batch, dtype=torch.float32)
        y = torch.tensor(y_batch, dtype=torch.float32)

        output = model(x)
        
        loss = criterion(output, y)
        loss.backward()
        optimizer.step()

        losses.append(float(loss))

        print(f"Epoch {i}: x={x}, y={y}, output={output}, loss={loss}")

    return np.array(losses)


def train_gsco2():
    
    model = models.gsCO2_model()
    loss = torch.nn.MSELoss()
    #opt = torch.optim.Adam(model.parameters(), lr=1e-2)
    opt = torch.optim.SGD(model.parameters(), lr=1e-2)
    epochs = 1000
    #data_iterator = gsco2_dummy_data_generator
    data_iterator = simple_data_generator

    losses = train_general(model, loss, opt, epochs, data_iterator)

    avg_window = 50
    avg_loss = losses[-avg_window:].mean()
    print(f"Average of the last {avg_window} losses: {avg_loss}")

    return model


"""
Represents an iterator objects that, for some function y = f(x),
samples random x values and outputs x,y pairs.
"""
class function_sample_iterator:  
    def __init__(self, f, sampler):
        self.f = f
        self.sampler = sampler
        
    def __iter__(self):
        return self

    def __next__(self):
        x = self.sampler()
        y = self.f(x)
        return x, y

gsco2_dummy_data_generator = function_sample_iterator(
    f = lambda x : np.array([x[0]*x[1]/((x[2]-x[3])*(1+x[4]/x[5]))]),
    sampler = lambda: np.random.normal(10, 2, 6) # NOTE: This could easily be changed
)

simple_data_generator = function_sample_iterator(
    f = lambda x : np.array([x[0] + 2*x[1] + 3*x[2] + 4*x[3] + 5*x[4] + 6*x[5]]),
    sampler = lambda: np.random.normal(0, 1, 6) # data centered around 0
)




