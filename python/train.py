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

import time

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

        # Obtain model predictions
        output = model(x)

        # Compute the loss
        loss = criterion(output, y)

        if torch.isnan(loss) or torch.isinf(loss) or (output==0.).any(): # Safety measure to avoid breaking training, shouldn't happen in practice
            if printIter:
                print(f"Found nan or inf loss or 0 output: don't compute gradient")
                print(f"loss={loss}, output={output}, y={y}")
            continue

        losses.append(float(loss))
        if printIter:
            loss_average = np.array(losses[-print_every_iter:]).mean()
            print(f"Epoch {i}: loss={loss}, average of last {print_every_iter} losses: {loss_average}")

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

    # split the data to avoid tuning with the final testing data
    data, _ = train_test_split(data, test_size=config.train.test_split, random_state=config.train.random_state)

    model_func = models.rs_model if config.pipeline.predict_rs else models.gsCO2_model
    sub_model_producer = lambda hps : model_func(hidden_size=hps['hidden_size'], n_hidden=hps['n_hidden'], 
                                        activation=config.model.activation, batch_norm=hps['batch_norm'])

    # choose a loss
    loss_criterion = torch.nn.L1Loss()

    # choose an optimizer
    opt_producer = lambda parameters, hps : torch.optim.Adam(parameters, lr=hps['lr'], weight_decay=hps['weight_decay'])

    # choose a data iterator
    iterator_producer = lambda data : batch_ctx_dict_iterator(data, batch_size=config.train.batch_size)

    # build the pipeline around the specific models
    def model_producer(hps):
        rs_model = sub_model_producer(hps)
        if config.pipeline.predict_rs:
            pipeline = pipelines.make_pipeline(rs_model, None, None, output_rs=False)
        else:
            pipeline = pipelines.make_pipeline(None, rs_model, None, output_rs=False)
        return pipeline, rs_model.parameters(), lambda: rs_model.eval()

    # perform K-Fold for the given set of hp
    # dict hp_name -> list of values
    param_possible_values = {'lr' : config.train.lr, 'weight_decay' : config.train.weight_decay, 
        'n_hidden' : config.model.n_hidden, 'hidden_size': config.model.hidden_size, 'batch_norm': config.model.batch_norm}

    # iterator that gives all possible combinations
    grid_search_hp_iterator = utils.grid_search_iterator(param_possible_values)

    lowest_score = np.inf
    best_hps = None
    losses = {}
    benchmarks = {}
    coeffs_determination = {}
    for hp_set in grid_search_hp_iterator:
        total_score, loss, benchmark, coeff_determination = k_folds(config.train, np.array(data), model_producer, opt_producer, loss_criterion, iterator_producer, hp_set)
        print(f"K-Fold total score for hp-set {hp_set}: {total_score}")
        if total_score < lowest_score:
            lowest_score = total_score
            best_hps = hp_set

        hp_set_string = utils.hp_set_to_string(hp_set)
        losses[hp_set_string] = loss
        benchmarks[hp_set_string] = benchmark
        coeffs_determination[hp_set_string] = coeff_determination

    print(f"Finished hp-tuning: found best hp-set = {best_hps}, yielding score {lowest_score}")

    # losses is a dict hp_set_string -> array of losses
    # benchmarks is a dict hp_set_string -> execution time 

    return best_hps, losses, benchmarks, coeffs_determination

def train_and_evaluate_pipeline(config, data):

    train_data, test_data = train_test_split(data, test_size=config.train.test_split, random_state=config.train.random_state)

    model = models.rs_model_from_config(config.model) if config.pipeline.predict_rs else models.gsCO2_model_from_config(config.model)

    # choose a loss
    #loss = torch.nn.MSELoss()
    loss_criterion = torch.nn.L1Loss()

    # choose an optimizer
    #opt = torch.optim.Adam(gsCO2_model.parameters() + Vmax_model.parameters(), lr=3e-4)
    #opt = torch.optim.SGD(gsCO2_model.parameters(), lr=lr, , weight_decay=weight_decay)
    opt = torch.optim.Adam(model.parameters(), lr=config.train.lr, weight_decay=config.train.weight_decay)

    # choose a data iterator
    # data_iterator = iter(train_data)
    #data_iterator = random_sample_iterator(train_data)
    data_iterator = batch_ctx_dict_iterator(train_data, batch_size=config.train.batch_size)

    # build the pipeline around the specific models
    if config.pipeline.predict_rs:
        model_wrapper = pipelines.make_pipeline(model, None, None, output_rs=False)
    else:
        model_wrapper = pipelines.make_pipeline(None, model, None, output_rs=False)

    # eval the model before training
    model.eval()
    eval_before_training = eval.eval_general(model_wrapper, train_data, loss_criterion)
    print(f"Eval before training: {eval_before_training}")

    # Init wandb
    run = None
    if config.wandb.use_wandb:
        run = wandb.init(
            project=config.wandb.project,
            tags=config.wandb.tags,
            config=OmegaConf.to_container(config, resolve=True),
        )

    # train the whole pipeline
    model.train()
    start = time.time()
    losses = train_general(model_wrapper, loss_criterion, opt, config.train.epochs, data_iterator, wandb_run=run)
    end = time.time()
    benchmark = end - start # training time in seconds

    # eval the model after training
    model.eval()
    eval_after_training = eval.eval_general(model_wrapper, test_data, loss_criterion)
    print(f"Eval after training: {eval_after_training}")

    eval_after_training_train = eval.eval_general(model_wrapper, train_data, loss_criterion)
    print(f"Eval after training on train data (to check overfitting): {eval_after_training_train}")

    # eval the empirical model (for comparaison)
    empirical_model_wrapper = pipelines.make_pipeline(None, None, None, output_rs=False)
    eval_empirical_model = eval.eval_general(empirical_model_wrapper, test_data, loss_criterion)
    print(f"Eval empirical model: {eval_empirical_model}")

    # Eval the R² score (for information)
    model_coeff = eval.coefficient_of_determination(model_wrapper, test_data)
    empirical_coeff = eval.coefficient_of_determination(empirical_model_wrapper, test_data)
    print(f"Coefficient of determination after training: {model_coeff}")
    print(f"Coefficient of determination of empirical model: {empirical_coeff}")

    # show info about loss (optional)
    plot.plot_losses([losses], labels=['Loss'], coefficients=[model_coeff], minmax=False, smoothing=0.1)

    # return rs_model, Vmax_model
    return model, losses, benchmark

# took inspiration from https://github.com/christianversloot/machine-learning-articles/blob/main/how-to-use-k-fold-cross-validation-with-pytorch.md
def k_folds(config, data, model_producer, opt_producer, criterion, iterator_producer, hp_set):

    k = config.k_folds
    k_fold_epochs = config.k_folds_epochs
    total_score = 0.
    kfold = KFold(n_splits=k, shuffle=True)

    losses = torch.tensor(0.)
    benchmarks = 0
    coeff_determination = 0

    for train_ids, test_ids in kfold.split(data):
        # then they use handy batch methods that I could use elsewhere also

        split_train_data = data[train_ids]
        iterator =  iterator_producer(split_train_data)
        model, parameters, eval_lambda = model_producer(hp_set)
        opt = opt_producer(parameters, hp_set)
        start = time.time()
        loss = train_general(model, criterion, opt, k_fold_epochs, iterator)
        end = time.time()
        losses += torch.tensor(losses)
        benchmarks += (end - start)
        with torch.no_grad():
            eval_lambda() # put the model in eval mode
            test_data_iterator = data[test_ids]
            # evaluate with our eval function the test data
            eval_score = eval.eval_general(model, test_data_iterator, criterion)
            print(f"Finished a fold, got result={eval_score}")
            coeff_determination += eval.coefficient_of_determination(model, test_data_iterator)
            total_score += eval_score
    
    # average the results over trains?
    return total_score / k, loss / k, benchmarks / k, coeff_determination / k


    
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
    batch_indices = np.random.choice(len(dict_array), batch_size, replace=False)
    batch = dict_array[batch_indices]

    return utils.make_batch(batch)

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
def train_rs(config):
    
    model = models.rs_model_from_config(config.model)
    def wrapper(x):
        x = x * torch.tensor([1e2, 1e-4, 1e-1, 1e1, 1e-1, 1e-2, 1]) # Optional: makes the wrapper take the actual predictor values as inputs
        model_output = model(x)
        rs_small = model_output
        rs = rs_small*1e2
        return rs
    loss = torch.nn.L1Loss()
    opt = torch.optim.Adam(model.parameters(), lr=0.005)
    epochs = 10000
    batch_size = 1
    data_iterator = batch_decorator_iterator(rs_dummy_data_generator, batch_size)

    losses = train_general(wrapper, loss, opt, epochs, data_iterator)

    avg_window = 50
    avg_loss = losses[-avg_window:].mean()
    print(f"Average of the last {avg_window} losses: {avg_loss}")

    plot.plot_losses([losses], ['dummy training loss'], minmax=False, smoothing=0.5)

    return wrapper

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
    
class batch_decorator_iterator:  
    def __init__(self, iterator, batch_size=32):
        self.it = iterator
        self.bs = batch_size
        
    def __iter__(self):
        return self

    # assumes x and y are tensors of shape m and 1.
    # groups them into shape (batch_size, m) and (batch_size, 1)
    def __next__(self):
        xs = []
        ys = []
        for _ in range(self.bs):
            x, y = next(self.it)
            xs.append(np.array(x))
            ys.append(np.array(y))

        xs = np.array(xs)
        ys = np.array(ys)
        return torch.tensor(xs), torch.tensor(ys)


def generate_points(sampler, count = 100, dtype=torch.float64):
    return [sampler() for _ in range(count)]

def empirical_model(x):
    a1=5
    go=10000
    An=x[0]
    Pre=x[1]
    Cc=x[2]
    GAM=x[3]
    Ds=x[4]
    Do=x[5]
    Ts=x[6]
    gsCO2 = go + a1*An*Pre/((Cc-GAM)*(1+Ds/Do)) ###  [umolCO2 / s m^2] -- Stomatal Conductance
    rsCO2=1/gsCO2 ### [ s m^2 / umolCO2 ] Stomatal resistence or Canopy 
    rsH20 = (rsCO2/1.64)*(1e6) ### [ s m^2 / molH20 ] Stomatal resistence or canopy 
    rs = rsH20*(273.15*Pre)/(0.0224*(Ts+273.15)*101325) ## [s/m]  Stomatal resistence or Canopy [convert stomatal resistence in terms of water volumes into another unit]
    return rs

def rs_sampler():
    lows =  np.array([2, 90000,  25, 1, 10  , 50, 5])
    highs = np.array([5, 95000, 50, 10,  300, 100, 20]) 

    return lambda: torch.tensor(np.random.uniform(lows, highs), dtype=torch.float32)

rs_dummy_data_generator = function_sample_iterator(
    f = lambda x : empirical_model(x),
    sampler = rs_sampler()
)

