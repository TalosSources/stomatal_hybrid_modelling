import torch
import numpy as np
import random # NOTE: could use a random method from np

from omegaconf import OmegaConf
import wandb

from sklearn.model_selection import KFold, train_test_split

import pipelines

import models
import plot
import eval
import utils

"""
General torch training method.
Args:
    model : callable with the x produced by the data_iterator
    criterion : callable with the y produced by data_iterator and the model output
    optimizer : torch optimizer containing model's parameters
    epochs : integer containing the training step count
    data_iterator : iterator yielding pairs (x, y)
    batch_size=1 : size of mini-batches
    print_every_iter : epoch count between two information prints during training
"""
def train_general(
    model,
    criterion,
    optimizer,
    epochs,
    data_iterator,
    print_every_iter=20,
    wandb_run=None
):

    losses = []

    for i in range(epochs):

        printIter = (i%print_every_iter) == 0

        # Reset the gradients
        optimizer.zero_grad()
        
        # produce a minibatch with x and y
        x, y = next(data_iterator)

        # TODO: Decision. In the paper, they ignore points where Q_LE < 0. Here, I clamp them to 0. A sensible choice could be to select valid points in batch making or even data loading (for the latter, infrastructure exists already)
        #if (y < 0).any():
        #    print(f"y == 0: skipping")
        #    continue
        #y = y.clamp(min=0)

        # Obtain model predictions
        output = model(x)

        # Compute the loss
        loss = criterion(output, y)

        if torch.isnan(loss) or torch.isinf(loss) or (output==0.).any(): # TODO: very crude fix. Actually understand why having output=0 breaks the whole pipeline. Also filter out and sanitize data properly
            if printIter:
                print(f"Found nan or inf loss or 0 output: don't compute gradient")
                print(f"loss={loss}, output={output}, y={y}")
            continue

        losses.append(float(loss))
        if printIter:
            loss_average = np.array(losses[-print_every_iter:]).mean()
            print(f"Epoch {i}: y={y}, output={output}, loss={loss}, average of last {print_every_iter} losses: {loss_average}")

        # compute the loss gradient
        loss.backward()

        # perform an optimization step
        optimizer.step()

        # log to wandb TODO: log other stuff: grad, ??? some stuff specific to this problem
        if wandb_run is not None:
            wandb_run.log({'loss' : float(loss)}, commit=True)

    return np.array(losses)

"""
Performs a search over the hyperparameter-space.
For each hp-set, test its performances with k-fold cross-validation
Return the best performing hp-set
"""
def perform_hp_tuning(config, data):

    # TODO: Since a lot of this is common to both methods, need to factor it

    rs_model_producer = lambda : models.rs_model(config.model)
    Vmax_model = None # NOTE: only train gsCO2 for now TO CONFIG whether to use it?

    # choose a loss TO CONFIG ? maybe, not necessary
    loss_criterion = torch.nn.L1Loss()

    # choose an optimizer TO CONFIG ? maybe, not necessary
    opt_producer = lambda parameters, hps : torch.optim.Adam(parameters, lr=hps['lr'], weight_decay=hps['weight_decay'])

    # choose a data iterator TO CONFIG ? maybe, not necessary
    iterator_producer = lambda data : batch_ctx_dict_iterator(data, batch_size=config.train.batch_size)

    # build the pipeline around the specific models
    def model_producer():
        rs_model = rs_model_producer()
        pipeline = pipelines.make_pipeline(rs_model, Vmax_model, output_rs=False)
        return pipeline, rs_model.parameters()

    # perform K-Fold for the given set of hp
    # For now, use the entire data to perform k-fold
    
    # dict hp_name -> list of values
    param_possible_values = {'lr' : config.train.lr, 'weight_decay' : config.train.weight_decay}

    # iterator that gives all possible combinations
    grid_search_hp_iterator = utils.grid_search_iterator(param_possible_values)

    lowest_score = np.inf
    best_hps = None
    for hp_set in grid_search_hp_iterator:
        total_score = k_folds(data, model_producer, opt_producer, loss_criterion, iterator_producer, hp_set)
        print(f"K-Fold total score = {total_score}")
        if total_score < lowest_score:
            lowest_score = total_score
            best_hps = hp_set

    print(f"Finished hp-tuning: found best hp = {best_hps}, yielding score {lowest_score}")

    # eval the model after training
    # rs_model.eval()
    # eval_after_training = eval.eval_general(model_wrapper, test_data[::eval_stride], loss_criterion)
    # print(f"Eval after training: {eval_after_training}")

    # eval the empirical model (for comparaison)
    # empirical_model_wrapper = pipelines.make_pipeline(None, None, output_rs=False)
    # eval_empirical_model = eval.eval_general(empirical_model_wrapper, test_data[::eval_stride], loss_criterion)
    # print(f"Eval empirical model: {eval_empirical_model}")

    # show info about loss 
    # TODO: Show info about gradient norm.
    #plot.plot_losses(losses)

    # return rs_model, Vmax_model
    return best_hps

def train_and_evaluate_pipeline(config, data):

    train_data, test_data = train_test_split(data, test_size=0.2) # TO CONFIG

    rs_model = models.rs_model(config.model)
    Vmax_model = None # NOTE: only train gsCO2 for now TO CONFIG whether to use it?

    # choose a loss TO CONFIG ? maybe, not necessary
    #loss = torch.nn.MSELoss()
    loss_criterion = torch.nn.L1Loss()

    # choose an optimizer TO CONFIG ? maybe, not necessary
    #opt = torch.optim.Adam(gsCO2_model.parameters() + Vmax_model.parameters(), lr=3e-4)
    #opt = torch.optim.SGD(gsCO2_model.parameters(), lr=lr, , weight_decay=weight_decay)
    opt = torch.optim.Adam(rs_model.parameters(), lr=config.train.lr, weight_decay=config.train.weight_decay)

    # choose a data iterator TO CONFIG ? maybe, not necessary
    # data_iterator = iter(train_data)
    #data_iterator = random_sample_iterator(train_data)
    data_iterator = batch_ctx_dict_iterator(train_data, batch_size=config.train.batch_size) # TODO: Make a batch iterator for the new pipeline

    # build the pipeline around the specific models
    model_wrapper = pipelines.make_pipeline(rs_model, Vmax_model, output_rs=False)

    eval_stride = 20 # NOTE: To config?
    # eval the model before training
    #rs_model.eval() # TODO: Would be better for the wrapper to offer this method. perhaps the wrapper should simply inherit from nn.module()
    #eval_before_training = eval.eval_general(model_wrapper, train_data, loss_criterion)
    #print(f"Eval before training: {eval_before_training}")

    # Init wandb
    # NOTE: Consider moving it in the other method, train_pipeline
    run = None
    if config.wandb.use_wandb:
        run = wandb.init(
            project=config.wandb.project,
            tags=config.wandb.tags,
            config=OmegaConf.to_container(config, resolve=True),
        )

    # train the whole pipeline
    rs_model.train()
    losses = train_general(model_wrapper, loss_criterion, opt, config.train.epochs, data_iterator, wandb_run=run)

    # eval the model after training
    rs_model.eval()
    eval_after_training = eval.eval_general(model_wrapper, test_data[::eval_stride], loss_criterion)
    print(f"Eval after training: {eval_after_training}")

    # eval the empirical model (for comparaison)
    empirical_model_wrapper = pipelines.make_pipeline(None, None, output_rs=False)
    eval_empirical_model = eval.eval_general(empirical_model_wrapper, test_data[::eval_stride], loss_criterion)
    print(f"Eval empirical model: {eval_empirical_model}")

    # show info about loss 
    # TODO: Show info about gradient norm.
    plot.plot_losses(losses)

    # return rs_model, Vmax_model
    return rs_model, None

# took inspiration from https://github.com/christianversloot/machine-learning-articles/blob/main/how-to-use-k-fold-cross-validation-with-pytorch.md
def k_folds(data, model_producer, opt_producer, criterion, iterator_producer, hp_set):

    k = 5 # TO CONFIG
    k_fold_epochs = 500 # TO CONFIG
    total_score = 0.
    kfold = KFold(n_splits=k, shuffle=True)

    for train_ids, test_ids in kfold.split(data):
        # then they use handy batch methods that I could use elsewhere also

        # simply train the pipeline using data[train_ids] or something NOTE: Need to reset stuff like model, optimizer, etc.
        split_train_data = data[train_ids]
        iterator =  iterator_producer(split_train_data)
        model, parameters = model_producer()
        opt = opt_producer(parameters, hp_set)
        train_general(model, criterion, opt, k_fold_epochs, iterator)
        with torch.no_grad():
            test_data_iterator = data[test_ids]
            # evaluate with our eval function the test data
            eval_score = eval.eval_general(model, test_data_iterator, criterion)
            print(f"Finished a fold, got result={eval_score}")
            total_score += eval_score
    
    # average the results over trains?
    return total_score / k


    
class batch_ctx_dict_iterator:
    def __init__(self, data, batch_size=16):
        # data is an array of dicts of params-values pairs
        grouped_data, probs = utils.group_by_ctx_id(data) # data is a dict ctx_id -> array of data points with the same context set
        self.grouped_data = grouped_data
        self.probs = probs
        self.n_groups = len(probs)
        self.batch_size = batch_size

    def __iter__(self):
        return self
    
    def __next__(self):
        # sample randomly a ctx_set_id
        id = np.random.choice(self.n_groups, p=self.probs)

        # then, sample a random batch with this specific ctx_set_id
        return ctx_dict_batch(self.grouped_data[id], self.batch_size)
      
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
    #batch = random.sample(dict_array, batch_size)
    batch_indices = np.random.choice(len(dict_array), batch_size, replace=False)
    batch = dict_array[batch_indices]

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


# ------------------------OLD-METHODS-USEFUL-FOR-THE-REPORT------------------------
def train_gsco2():
    
    model = models.rs_model()
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
