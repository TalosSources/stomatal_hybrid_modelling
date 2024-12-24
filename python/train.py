import torch
import numpy as np

import photosynthesis_biochemical
import pipelines

import models
import plot
import utils
import eval

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

        printIter = (i%50) == 0
        
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

        # TODO: Decision. In the paper, they ignore points where Q_LE < 0. Here, I clamp them to 0.
        #if (y < 0).any():
        #    print(f"y == 0: skipping")
        #    continue
        #y = y.clamp(min=0)

        output, rs_values = model(x)
        #print(f"in train: output = {output}")
        #output = output[0]
        loss = criterion(output, y)
        #print(f"Epoch {i}: x={x}, y={y}, output={output}, loss={loss}")
        if printIter:
            print(f"Epoch {i}: y={y}, output={output}, loss={loss}")
        #print(f"Epoch {i}: loss={loss}")
        if torch.isnan(loss) or torch.isinf(loss) or (output==0.).any(): # TODO: very crude fix. Actually understand why having output=0 breaks the whole pipeline. Also filter out and sanitize data properly
            if printIter:
                print(f"Found nan or inf loss: don't compute gradient")
            continue
            #return

        # retain grads
        for rs in rs_values:
            rs.retain_grad()
        loss.retain_grad()
        output.retain_grad()
        for p in optimizer.param_groups:
            for param in p['params']:
                #print(f"retaining grad on param {param} (and data={param.data})")
                param.retain_grad()

        # backward
        loss.backward()

        # print infos
        #utils.printGradInfo(output, "outputAfterBackward")
        #for rs in rs_values:
        #    utils.printGradInfo(rs, "rsInTrainGeneral")
        #utils.printGradInfo(loss, "loss")
        #utils.printGradInfo(photosynthesis_biochemical.gsCO2, "gsco2")
        #utils.printGradInfo(photosynthesis_biochemical.model_output, "gsco2ModelOutput")
        #for p in optimizer.param_groups:
        #    for param in p['params']:
        #        utils.printGradInfo(param, "paramDataX")

        # step
        optimizer.step()

        losses.append(float(loss))

        #print(f"----------------------------------------------------")

        

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
    lr = 1e-2
    epochs = 10000

    model_wrapper = pipelines.make_pipeline(gsCO2_model, Vmax_model, output_rs=True)
    data_iterator = batch_ctx_dict_iterator(train_data, batch_size=32) # TODO: Make a batch iterator for the new pipeline
    # data_iterator = iter(train_data)
    #data_iterator = random_sample_iterator(train_data)

    # choose a loss
    #loss = torch.nn.MSELoss()
    loss = torch.nn.L1Loss()

    # choose an optimizer
    #opt = torch.optim.Adam(gsCO2_model.parameters() + Vmax_model.parameters(), lr=3e-4)
    opt = torch.optim.Adam(gsCO2_model.parameters(), lr=lr)
    #opt = torch.optim.SGD(gsCO2_model.parameters(), lr=lr)

    print(f"using parameters: {gsCO2_model.parameters()}")
    #parameters = gsCO2_model.parameters()
    print(f"Gsco2 Model require_grad: {next(gsCO2_model.parameters()).requires_grad}")
    for p in gsCO2_model.parameters():
        print(f"param p : {p}")
        print(f"requires_grad? = {p.requires_grad}")

    # eval the model before training
    eval_before_training = eval.eval_general(model_wrapper, train_data, loss)
    print(f"Eval before training: {eval_before_training}")

    losses = train_general(model_wrapper, loss, opt, epochs, data_iterator)

    # eval the model after training
    eval_after_training = eval.eval_general(model_wrapper, train_data, loss)
    print(f"Eval after training: {eval_after_training}")

    # eval the empirical model
    empirical_model_wrapper = pipelines.make_pipeline(None, None, output_rs=True)
    eval_empirical_model = eval.eval_general(empirical_model_wrapper, train_data, loss)
    print(f"Eval empirical model: {eval_empirical_model}")

    # show info about loss NOTE: Consider using wandb or similar foss
    # TODO: Show info about gradient norm.
    avg_window = 50
    avg_loss = losses[-avg_window:].mean()
    print(f"Average of the last {avg_window} losses: {avg_loss}")
    plot.plot_losses(losses)

    return gsCO2_model, Vmax_model




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
    
class batch_ctx_dict_iterator:
    def __init__(self, data, batch_size=16):
        # data is an array of dicts of params-values pairs
        self.data = data
        self.batch_size = batch_size

    def __iter__(self):
        return self
    
    def __next__(self):
        return ctx_dict_batch(self.data, self.batch_size)
      
class random_sample_iterator:
    def __init__(self, data):
        # data is an array of (input, output) tuples, where input and output can be anything.
        self.data = data
        self.n = len(data)

    def __iter__(self):
        return self
    
    def __next__(self):
        idx = random.randint(0, self.n-1)
        #print(f"(just sampled {idx})")
        return self.data[idx]

    
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

"""
This function takes a large array where each line is a tuple (x, y), 
where x is a tuple (ctx_dict, global_dict) and y is a tensor with a single value,
where ctx_dict is a dict ctx -> predictor_dict,
and where predictor_dict and global_dict are dicts key -> tensor with a single value.
and return a random batch of batch_size, in the form of a tuple (x, y) where 
y is an tensor of size batch_size, and x has the same form, except
predictor_dict and global_dict map keys to tensors of size batch_size
Of course, the elements must be in the same order in the tensor for each key.
"""
def ctx_dict_batch(dict_array, batch_size):

    # batch is an array of size batch_size containing (x,y) tuples
    batch = random.sample(dict_array, batch_size)

    # predictors is an array of size batch_size of (ctx_dict, global_dict) values
    # y is an array of size batch_size of y values
    predictors, outputs = zip(*batch)

    # outputs_batch is already a y tensor in the form we want it
    outputs_batch = torch.tensor(outputs) 

    # ctx_dicts is an array of size batch_size of ctx_dict values, same for global_dicts
    ctx_dicts, global_dicts = zip(*predictors)

    # for each key present in global_dicts, we collect all the batch_size values for this key in a single tensor
    # and map the key to that tensor
    global_batch = {key: torch.tensor([global_dict[key] for global_dict in global_dicts]) 
                  for key in global_dicts[0].keys()}
    
    # for each ctx present in ctx_dicts, and for each key present for this ctx,
    # aggregate all these values into a tensor, and map the corresponding key to this tensor,
    # and the corresponding ctx to this predictor_dict. 
    ctx_batch = {
        ctx: {
            key: torch.tensor([ctx_dict[ctx][key] for ctx_dict in ctx_dicts]) 
            for key in ctx_dicts[0][ctx].keys()
        }
        for ctx in ctx_dicts[0].keys()
    }

    # ugly fix TODO I'm not sure what the proper way of doing this would be 
    for predictor_dict in ctx_batch.values():
        predictor_dict['CT'] = 3

    # TODO: Sort out the reason that output and y don't have the same dimension

                  
    return ((ctx_batch, global_batch), outputs_batch)

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




