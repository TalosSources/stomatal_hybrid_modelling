from photosynthesis_biochemical import photosynthesis_biochemical
import differentiable_relations

import models
import torch
import numpy as np

import plot
import inspect

def train_general(
    model,
    criterion,
    optimizer,
    epochs,
    data_iterator,
    batch_size=1
):

    losses = []

    for i in range(epochs):
        
        optimizer.zero_grad()
        

        if batch_size > 1:
            x_batch = []
            y_batch = []
            for _ in range(batch_size):
                x, y = next(data_iterator)
                x_batch.append(x)
                y_batch.append(y)
            x = torch.tensor(np.array(x_batch), dtype=torch.float32)
            y = torch.tensor(np.array(y_batch), dtype=torch.float32)
        else:
            x, y = next(data_iterator)

        output = model(x)
        loss = criterion(output, y)
        #print(f"Epoch {i}: x={x}, y={y}, output={output}, loss={loss}")
        print(f"Epoch {i}: y={y}, output={output}, loss={loss}")
        #print(f"Epoch {i}: loss={loss}")
        if torch.isnan(loss) or torch.isinf(loss) or (output==0.).any(): # TODO: very crude fix. Actually understand why having output=0 breaks the whole pipeline. Also filter out and sanitize data properly
            print(f"Found nan or inf loss: don't compute gradient")
            continue
            #return
        #loss.backward()
        #optimizer.step()

        losses.append(float(loss))

        

    return np.array(losses)


def train_gsco2():
    
    model = models.gsCO2_model()
    loss = torch.nn.MSELoss()
    opt = torch.optim.Adam(model.parameters(), lr=3e-4)
    #opt = torch.optim.SGD(model.parameters(), lr=1e-2)
    epochs = 10000
    data_iterator = gsco2_dummy_data_generator
    #data_iterator = simple_data_generator

    losses = train_general(model, loss, opt, epochs, data_iterator)

    avg_window = 50
    avg_loss = losses[-avg_window:].mean()
    print(f"Average of the last {avg_window} losses: {avg_loss}")

    plot.plot_losses(losses)
    #plot.plot_univariate_slice(model, gsco2_dummy_data_generator.f, torch.tensor([5, 5, 8, 2, 2, 4], dtype=torch.float32), 2, np.arange(5, 11, 0.1), "gsco2")
    #plot.plot_univariate_slice(model, gsco2_dummy_data_generator.f, torch.tensor([5, 5, 8, 2, 2, 4], dtype=torch.float32), 1, np.arange(2, 8, 0.1), "gsco2")
    #plot.plot_univariate_slice(model, gsco2_dummy_data_generator.f, torch.tensor([5, 5, 8, 2, 2, 4], dtype=torch.float32), 0, np.arange(2, 8, 0.1), "gsco2")
    plot.fit_plot(model, gsco2_dummy_data_generator.f, generate_points(gsco2_dummy_data_generator.sampler, 100), "gsco2_fcn")

    return model

def train_vm():
    
    model = models.vm_model()
    loss = torch.nn.MSELoss()
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    epochs = 500
    data_iterator = vm_dummy_data_generator

    losses = train_general(model, loss, opt, epochs, data_iterator, batch_size=1)

    avg_window = 50
    avg_loss = losses[-avg_window:].mean()
    print(f"Average of the last {avg_window} losses: {avg_loss}")

    plot.plot_losses(losses)
    plot.plot_univariate_slice(model, vm_dummy_data_generator.f, torch.tensor([5, 5], dtype=torch.float32), 0, np.arange(3, 7, 0.1), "vm")
    plot.plot_univariate_slice(model, vm_dummy_data_generator.f, torch.tensor([5, 5], dtype=torch.float32), 1, np.arange(3, 7, 0.1), "vm")
    plot.fit_plot(model, vm_dummy_data_generator.f, generate_points(vm_dummy_data_generator.sampler, 100), "vm_fcn")

    return model

def train_pipeline(train_data):
    
    gsCO2_model = models.gsCO2_model()
    #Vmax_model = models.vm_model()
    Vmax_model = None # NOTE: only train gsCO2 for now

    """
    Need some function that computes pb outputs given chosen predictors. 
    It should store the constant inputs not considered predictors. 
    Else we should create another train_general function 
    """
    model_wrapper = make_pipeline(gsCO2_model, Vmax_model)
    #data_iterator = batch_dict_iterator(train_data, batch_size=1)
    data_iterator = iter(train_data)

    loss = torch.nn.MSELoss()
    #opt = torch.optim.Adam(gsCO2_model.parameters() + Vmax_model.parameters(), lr=3e-4)
    opt = torch.optim.Adam(gsCO2_model.parameters(), lr=1e-2)

    epochs = 4000

    losses = train_general(model_wrapper, loss, opt, epochs, data_iterator)

    # show info about loss NOTE: Consider using wandb or similar foss
    avg_window = 50
    avg_loss = losses[-avg_window:].mean()
    print(f"Average of the last {avg_window} losses: {avg_loss}")
    plot.plot_losses(losses[100:])

    return gsCO2_model, Vmax_model

# TODO: Move pure pipeline stuff (i.e. mostly agnostic to training) in a separate pipeline file
def make_pipeline(gsCO2_model, Vmax_model, output_rs=False):

    # define constants. Doing our function structure that way allows to store the constants cleanly
    # constants = ... # NOTE: Might not be necessary

    # the pipeline calls the pb module with the 2 models inserted, using the constants and predictors as inputs, and converts the outputs to Q_LE 
    pb_params = set(inspect.signature(photosynthesis_biochemical).parameters)
    qle_params = set(inspect.signature(differentiable_relations.Q_LE).parameters)

    def subpipeline(predictors):
        pb_predictors = {k: v for k, v in predictors.items() if k in pb_params}
        _,_,rs,_,_,_,_ = photosynthesis_biochemical(**pb_predictors, gsCO2_model=gsCO2_model, Vmax_model=Vmax_model)
        #print(f"input_rs={input_rs} while output_rs = {rs}")
        qle_predictors = {k: v for k, v in predictors.items() if k in qle_params}
        qle_predictors.__delitem__("rs") # TODO: temporary debug fix

        pred_rs = predictors["rs"]
        print(f"given_rs: {pred_rs}, computed_rs: {rs}")
        if not torch.isnan(pred_rs) and not torch.isinf(pred_rs) and not torch.isnan(rs) and not torch.isinf(rs):
            if torch.abs(rs - pred_rs) > 0.5:
                print(f"DIFFERENT!")
                print(f"and pb_predictors are: {pb_predictors}")

        # COMPUTE Q_LE USING THE PREDICTOR RS
        #Q_LE = differentiable_relations.Q_LE(rs=pred_rs, **qle_predictors)

        # COMPUTE Q_LE USING THE ESTIMATED RS
        Q_LE = differentiable_relations.Q_LE(rs=rs, **qle_predictors)

        #return Q_LE, rs if output_rs else Q_LE
        #print(f"\n\nsubpipeline got {predictors}\n\n and output {Q_LE} (shape={Q_LE.size()})\n\n")
        return Q_LE
    
    def W_on_m2_to_mm_on_H(Q):
        # For now simply divide by 684 NOTE: Not 100% sure it's the way to go
        return Q / 684
    
    def mm_on_H_to_W_on_m2(Q):
        return Q * 684
    
    """
    This is the differentiable pipeline, taking as input an object containing
    values for the 4 different contexts, computing the Q_LE term for each of 
    them, and aggregating them to obtain a Q_LE comparable with observations
    Assumes it receives predictors p, on this format:
    * single value tensors for general EG, ELitter, etc., and also EIn_H and EIn_L 
    * For key ctx, a map of predictors for pb and Q_LE, with values for this ctx.
    """
    def pipeline(p):
        single_pipeline_values, g, output = p
        # predictors contains a map between context and predictor maps for this context.
        Q_LE_wm2 = torch.zeros(1)
        #print(f"the acc is {Q_LE_wm2}, size={Q_LE_wm2.size()}")
        for ctx, preds in single_pipeline_values.items():
            print(f"[ctx={ctx}] ", end='')
            ctx_res = subpipeline(preds)
            if not torch.isnan(ctx_res):
                #print(f"... adding it ...")
                Q_LE_wm2 += ctx_res
        Q_LE_mmH = W_on_m2_to_mm_on_H(Q_LE_wm2)
        #print(f"predicted Q_LE_mmH: {Q_LE_mmH}")

        HLTotal = Q_LE_mmH + g["EIn_H"] + g["EIn_L"]
        ET = HLTotal + g["EG"] + g["ELitter"] + g["ESN"] + g["ESN_In"] + g["EWAT"] + g["EICE"] + g["EIn_urb"] + g["EIn_rock"]

        final_Q_LE = mm_on_H_to_W_on_m2(ET)

        return final_Q_LE

    return pipeline



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
    
class batch_dict_iterator:
    def __init__(self, data, batch_size=16):
        # data is an array of dicts of params-values pairs
        self.data = data
        self.batch_size = batch_size

    def __iter__(self):
        return self
    
    def __next__(self):
        return dict_batch(self.data, self.batch_size)

    
def generate_points(sampler, count = 100):
    return [torch.tensor(sampler(), dtype=torch.float32) for _ in range(count)]

gsco2_dummy_data_generator = function_sample_iterator(
    f = lambda x : np.array([x[0]*x[1]/((x[2]-x[3])*(1+x[4]/x[5]))]),
    sampler = lambda: np.random.normal([5, 5, 8, 2, 2, 4], 1, 6) # NOTE: This could easily be changed
)

def vm_empirical(x):
    # input in x
    Ts = torch.tensor(x[0], dtype=torch.float32)
    Vmax = torch.tensor(x[1], dtype=torch.float32)

    # constants?
    Ts_k = torch.tensor(Ts + 273.15) ##[K]
    Tref = torch.tensor(25 + 273.15) ## [K] Reference Temperature
    Ha = torch.tensor(65.33) ##  [kJ/mol] Activation Energy - Plant Dependent
    DS = torch.tensor(0.485) ## [kJ / mol K]  entropy factor - Plant Dependent
    R =   torch.tensor(0.008314) ##  [kJ�/ K mol] Gas Constant
    Hd =  torch.tensor(200) # [kJ/mol]  Deactivation Energy -- constant

    kT =torch.exp(Ha*(Ts_k-Tref)/(Tref*R*Ts_k))*(1+torch.exp((Tref*DS - Hd)/(Tref*R)))/(1+torch.exp((Ts_k*DS-Hd)/(Ts_k*R)))
    Vm=Vmax*kT

    return Vm
    #return Ts + Vmax

    # in some strange case
    #Vm = Vmax*fT*f1T*f2T

vm_dummy_data_generator = function_sample_iterator(
    f = vm_empirical,
    sampler= lambda: np.random.normal(5, 1, 2)
)

simple_data_generator = function_sample_iterator(
    f = lambda x : np.array([x[0] * x[1]]),# + x[2] + x[3] + x[4] + x[5]]),
    sampler = lambda: np.random.normal(0, 2, 2) # data centered around 0
)


"""
This function takes a large array where each line is a dict with values for parameter names, 
and return a random batch of batch_size, in the form of a dict where the value for a 
parameter name is a tensor of size batch_size.
Of course, the elements must be in the same order in the tensor for each parameter.
"""
import random
def dict_batch(dict_array, batch_size):
    batch = random.sample(dict_array, batch_size)
    predictors, outputs = zip(*batch)
    predictors_batch = {key: torch.tensor([entry[key] for entry in predictors]) 
                  for key in predictors[0].keys()}

    # ugly fix TODO
    predictors_batch['CT'] = 3

    # TODO: Sort out the reason that output and y don't have the same dimension

    outputs_batch = torch.tensor(outputs) 
    #outputs_batch[:] = 7.58 # dummy output, just to see if the training system can learn this TODO: remove
                  
    return (predictors_batch, outputs_batch)


def check_gradients(module : torch.nn.Module):
    for name, param in module.named_parameters():
        if param.grad is not None:
            if torch.isnan(param.grad).any():
                print(f"NaN detected in gradients of {name}")
            elif torch.isinf(param.grad).any():
                print(f"Inf detected in gradients of {name}")
    for param in module.parameters():
        if param.grad is not None:
            if torch.isnan(param.grad).any():
                print(f"NaN detected in gradients of {param}")
            elif torch.isinf(param.grad).any():
                print(f"Inf detected in gradients of {param}")




