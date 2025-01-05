"""
This file contains method useful for evaluating, interpreting and plotting 
the results and components of the experiments.
This is used in writing the report for this project.
"""
import numpy as np
import torch
import os
from omegaconf import OmegaConf

from sklearn.model_selection import train_test_split

from copy import deepcopy

import torchviz

import plot
import pipelines
import models
import data
import differentiable_relations

import train

# GENERAL TODO: Be careful with colors. Try to use the same consistent colors for the same things, like rs_model, etc. 
# Try to use sensible colors, like red for hot, blue for cold  

def load_model(config_name):
    # load the model weights
    results_path = os.path.join('results', config_name)
    rs_model_state_dict = torch.load(os.path.join(results_path, "model_weights.pt"))
    cfg = OmegaConf.load(f"configs/{config_name}.yaml")
    rs_model = models.rs_model_from_config(cfg.model)
    rs_model.load_state_dict(rs_model_state_dict)
    return rs_model

def load_train_test_data(config_name):
    config = OmegaConf.load(f"configs/{config_name}.yaml")
    all_data = data.load_pipeline_data_dict_from_all_sites(config.data, verbose=True, output_keys=['LE']) # NOTE: LE_CORR?
    train_data, test_data = train_test_split(all_data, test_size=config.train.test_split, random_state=config.train.random_state)
    return train_data, test_data

def load_ordered_data(config_name):
    config = OmegaConf.load(f"configs/{config_name}.yaml")
    site_names = config.data.sites
    base_path = os.path.expanduser(config.data.base_path)
    combined_site_data = []
    for site_name in site_names:
        pred_path = os.path.join(base_path, f"Results_{site_name}.mat")
        obs_path = os.path.join(base_path, f"Res_{site_name}.mat")
        exclude_negative_outputs = False
        site_data = []
        n_points = 20000 #1000 corresponds to about 3 weeks NOTE: Where should we choose this?
        data.load_pipeline_data_dict_single_site(pred_path, obs_path, site_data, 
                                             nPoints=n_points, exclude_negative_outputs=exclude_negative_outputs, verbose=True)
        if len(site_data) > 0:
            combined_site_data.append((site_data, site_name))
    return combined_site_data
    
# TODO Re-run hp-tuning with the new configs, and choose the new interesting runs to plot
# TODO: Replace R^2 with the actual coefficients
# TODO: Latex table with coefficients? or perhaps another tool for the table
def generate_hp_tuning_loss_plot():
    results_path = "results/lr_wd_tuning"
    filenames = [
        #"losses_lr-0.01_weight_decay-0_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        #"losses_lr-0.001_weight_decay-0_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        #"losses_lr-0.0001_weight_decay-0_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        "losses_lr-0.01_weight_decay-0.0001_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        "losses_lr-0.001_weight_decay-0.0001_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        "losses_lr-0.0001_weight_decay-0.0001_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        #"losses_lr-0.01_weight_decay-1e-06_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        #"losses_lr-0.001_weight_decay-1e-06_n_hidden-4_hidden_size-128_batch_norm-True_.npy",
        #"losses_lr-0.0001_weight_decay-1e-06_n_hidden-4_hidden_size-128_batch_norm-True_.npy"
    ]

    losses = [
        np.load(os.path.join(results_path, filename)) for filename in filenames
    ]

    labels = [
        #"lr=1e-2, lambda=0",
        #"lr=1e-3, lambda=0",
        #"lr=1e-4, lambda=0",
        "lr=1e-2, lambda=1e-4",
        "lr=1e-3, lambda=1e-4",
        "lr=1e-4, lambda=1e-4",
        #"lr=1e-2, lambda=1e-6",
        #"lr=1e-3, lambda=1e-6",
        #"lr=1e-4, lambda=1e-6",
    ]

    coefficients = [
        0.25,
        0.36,
        0.78
    ]

    plot.plot_losses(losses, labels, coefficients=coefficients, minmax=False, smoothing=0.01)

# TODO: Choose interesting x_0s
# TODO: Make a common plot with all chosen variables? plt.subplot... requires factoring everything...
def generate_rs_slice_plots(rs_model):
    # load specific(s) datapoint(s)
    # as info: predictors = torch.stack([An * 1e2,Pre * 1e-4,Cc * 1e-1,GAM * 1e1,Ds * 1e-1,Do * 1e-2, Ts], dim=-1)
    # Given: Pre, Cc, Ds, Ts, Do. An and GAM are computed
    # It is very annoying and ugly to get sensible values for x_0 and sensible variation intervals in a clean way.
    # Either I copy paste the entire pb.py file, or output many variables in pb just for this...
    # The simplest way would be to print the predictors just before passing them. Then, I can vary the ones that make the most sense.
    # I'll vary Ds, Ts, Cc, and Pre.
    # An=[-0.2730], Pre=[91887.5000], Cc=[38.0901], GAM=[2.2750], Ds=[393.0414], Do=[1000.], Ts=[14.6980]
    with torch.no_grad():
        x_0s = [
            torch.tensor([-0.2730, 91887.5000, 38.0901, 2.2750, 393.0414, 1000., 14.6980]),
            #torch.tensor([2.510115,   84400,   24.540,  2.2750,  898.362,   800.00,  12.1075])
        ]

        x_0_labels = [
            '', 
            #''
        ]

        # wrap the rs_model so that it outputs the actual rs
        def wrapper(x):
            x = x * torch.tensor([1e2, 1e-4, 1e-1, 1e1, 1e-1, 1e-2, 1]) # Optional: makes the wrapper take the actual predictor values as inputs
            model_output = rs_model(x)
            rs_small = torch.nn.functional.softplus(model_output)
            return rs_small * 1e2
        
        def empirical_model(x):
            a1=5
            go=10000
            An=x[:,0]
            Pre=x[:,1]
            Cc=x[:,2]
            GAM=x[:,3]
            Ds=x[:,4]
            Do=x[:,5]
            Ts=x[:,6]
            gsCO2 = go + a1*An*Pre/((Cc-GAM)*(1+Ds/Do)) ###  [umolCO2 / s m^2] -- Stomatal Conductance
            gsCO2 = torch.max(gsCO2, torch.tensor(go))
            rsCO2=1/gsCO2 ### [ s m^2 / umolCO2 ] Stomatal resistence or Canopy 
            rsH20 = (rsCO2/1.64)*(1e6) ### [ s m^2 / molH20 ] Stomatal resistence or canopy 
            rs = rsH20*(273.15*Pre)/(0.0224*(Ts+273.15)*101325) ## [s/m]  Stomatal resistence or Canopy [convert stomatal resistence in terms of water volumes into another unit]
            return rs

        # Define what to plot
        idx_and_ranges = [
            [1, (82000, 102000), 'Pressure [Pa]'], # Pre plot
            [2, (25, 50), 'Cc [Pa*molCO2/molAIR]'], # Cc plot
            [4, (0, 4000), 'VPD [Pa]'], # Ds plot
            [6, (-10, 40), 'Temperature [°C]'], # Ts plot
        ]

        # we might add more models on the same plot
        labels = ['rs_model (ours)', 'empirical model']
        models_list = [wrapper, empirical_model]

        # for each of the input variables, make a slice plot
        for idx, idx_range, idx_label in idx_and_ranges:
            path = os.path.join('figures', 'rs_slices', f'slice_varying_{idx}.png')
            plot.plot_univariate_slices(models_list, x_0s, x_0_labels, idx, idx_label, idx_range, 100,
                                        labels, res_label='rs [s/m]',path=path, show=True)

# TODO: Choose interesting x_0s (and datapoint)
def generate_Q_LE_slice_plots(rs_model):
    # Idea: work with lambda. have a lambda closure that contains 'hardcoded' (using fetched data) predictors for all values except
    # the ones that we vary
    # load specific(s) datapoint(s)
    # I'll vary Ds, Ts, Cc, and Pre.
    # An=[-0.2730], Pre=[91887.5000], Cc=[38.0901], GAM=[2.2750], Ds=[393.0414], Do=[1000.], Ts=[14.6980]
    with torch.no_grad():

        t = torch.tensor
        datapoint = (
            (
                {
                    'sun_H': {'Cc': t([368.0680]), 'IPAR': t([0.]), 'Csl': t([362.6800]), 'ra': t([38.7362]), 'rb': t([19.6884]), 'Ts': t([-5.3135]), 'Pre': t([823.]), 'Ds': t([50.7378]), 'Psi_L': t([0.]), 'Rn': t([-58.5669]), 'QG': t([-68.8654]), 'Vmax': t([25.3091]), 'Psi_sto_50': t([-2.5000]), 'Psi_sto_00': t([-0.5000]), 'CT': t([3.]), 'Ha': t([72.]), 'FI': t([0.0810]), 'Oa': t([210000.]), 'Do': t([800.]), 'a1': t([5.]), 'go': t([0.0100]), 'gmes': t([np.inf]), 'rjv': t([2.1000]), 'DS': t([0.6490])},
                    #'sun_L': {'Cc': t([368.0680]), 'IPAR': t([0.]), 'Csl': t([362.6800]), 'ra': t([38.7362]), 'rb': t([19.6884]), 'Ts': t([-5.3135]), 'Pre': t([823.]), 'Ds': t([50.7378]), 'Psi_L': t([0.]), 'Rn': t([-58.5669]), 'QG': t([-68.8654]), 'Vmax': t([25.3091]), 'Psi_sto_50': t([-2.5000]), 'Psi_sto_00': t([-0.5000]), 'CT': t([3.]), 'Ha': t([72.]), 'FI': t([0.0810]), 'Oa': t([210000.]), 'Do': t([800.]), 'a1': t([5.]), 'go': t([0.0100]), 'gmes': t([np.inf]), 'rjv': t([2.1000]), 'DS': t([0.6490])},
                    #'shd_H': {'Cc': t([368.0680]), 'IPAR': t([0.]), 'Csl': t([362.6800]), 'ra': t([38.7362]), 'rb': t([19.6884]), 'Ts': t([-5.3135]), 'Pre': t([823.]), 'Ds': t([50.7378]), 'Psi_L': t([0.]), 'Rn': t([-58.5669]), 'QG': t([-68.8654]), 'Vmax': t([25.3091]), 'Psi_sto_50': t([-2.5000]), 'Psi_sto_00': t([-0.5000]), 'CT': t([3.]), 'Ha': t([72.]), 'FI': t([0.0810]), 'Oa': t([210000.]), 'Do': t([800.]), 'a1': t([5.]), 'go': t([0.0100]), 'gmes': t([np.inf]), 'rjv': t([2.1000]), 'DS': t([0.6490])},
                    #'shd_L': {'Cc': t([368.0680]), 'IPAR': t([0.]), 'Csl': t([362.6800]), 'ra': t([38.7362]), 'rb': t([19.6884]), 'Ts': t([-5.3135]), 'Pre': t([823.]), 'Ds': t([50.7378]), 'Psi_L': t([0.]), 'Rn': t([-58.5669]), 'QG': t([-68.8654]), 'Vmax': t([25.3091]), 'Psi_sto_50': t([-2.5000]), 'Psi_sto_00': t([-0.5000]), 'CT': t([3.]), 'Ha': t([72.]), 'FI': t([0.0810]), 'Oa': t([210000.]), 'Do': t([800.]), 'a1': t([5.]), 'go': t([0.0100]), 'gmes': t([np.inf]), 'rjv': t([2.1000]), 'DS': t([0.6490])},
                }, 
                {'EIn_H': t([0.]), 'EIn_L': t([0.]), 'EG': t([0.0028]), 'ELitter': t([0.]), 'ESN': t([0.]), 'ESN_In': t([0.]), 'EWAT': t([0.]), 'EICE': t([0.]), 'EIn_urb': t([0.]), 'EIn_rock': t([0.])}
            ), 
            t(2.6774)
        )

        x_0s = [torch.tensor([823. , 50.7378, -5.3135])]
        x_0_labels = ['dummy point']
        
        pipeline = pipelines.make_pipeline(rs_model, None, None)
        empirical_pipeline = pipelines.make_pipeline(None, None, None)

        steps = 100

        def wrapper(pip, x):
            # x is [Pre, Ds, Ts]
            datapoint_2 = deepcopy(datapoint)
            p, _ = datapoint_2
            ctx_preds, _ = p
            for _, pred in ctx_preds.items(): # NOTE: Does it vary datapoint_2?
                pred['Pre'] = x[:,0]
                pred['Ds'] = x[:,1]
                pred['Ts'] = x[:,2]
                for k, v in pred.items():
                    if k not in ['Pre', 'Ds', 'Ts']:
                        pred[k] = pred[k].expand(steps)
            return pip(p)

        our_wrapper = lambda x : wrapper(pipeline, x)
        empirical_wrapper = lambda x : wrapper(empirical_pipeline, x)

        idx_and_ranges = [
            [0, (820, 1020), 'Pressure [Pa]'], # Pre plot
            [1, (0, 4000), 'VPD [Pa]'], # Ds plot
            [2, (-10, 40), 'Temperature [°C]'], # Ts plot
        ]

        # we might add more models on the same plot
        labels = ['Q_LE pipeline (ours)', 'empirical model']
        models_list = [our_wrapper, empirical_wrapper]

        # for each of the input variables, make a slice plot
        for idx, idx_range, idx_label in idx_and_ranges:
            path = os.path.join('figures', 'Q_LE_slices', f'slice_varying_{idx}.png')
            plot.plot_univariate_slices(models_list, x_0s, x_0_labels, idx, idx_label, idx_range, steps, labels, res_label='Q_LE [W/m²]', path=path, show=True)

# TODO: Do daily aggregates. Idea: log in global_preds, the i index of the timestep. 
def generate_scatter_plots(test_data, rs_model):
    # then, write a function that, given data from some site (in sequential order? perhaps not needed),
    # groups the data in  
    
    test_data = test_data[:5000] # take only a subset of the data. for performances. TODO OPT: Use ctx batches in the plot method (or here?), for better perf because it's so slow now.

    # make the pipeline
    pipeline = pipelines.make_pipeline(rs_model, None, None)
    print(f"finished making pipeline")

    empirical_pipeline = pipelines.make_pipeline(None, None, None)

    models_list = [pipeline, empirical_pipeline]
    labels = ['our pipeline', 'empirical pipeline']

    path = os.path.join('figures', 'scatter_plots', f'simple_best_model_scattering.png')

    with torch.no_grad():
        plot.fit_plot(models_list, test_data, labels, path, plot_range=(0, 600), show=True)

# TODO: Put annotations on the plot I choose
def generate_torch_viz_graph():
    # First, decide if I want to keep the model architecture. Perhaps not, it's not the focus
    # Then, give some of the predictors to torchviz, to display the actual Q_LE PM pipeline
    # If possible, find a datapoint with 2 or 4 valid contexts. (further? perhaps the crazy AU-Stp?)
    # Then, try to remove that long line of ADDs and stuff.
    # Finally, try to give better names to the nodes if possible

    cfg = OmegaConf.load("configs/default.yaml")
    cfg.data.nPoints = 10
    cfg.data.sites = ["CH-Dav"]
    cfg.model.batch_norm = False
    rs_model = models.rs_model_from_config(cfg.model)

    class BlackBox(torch.nn.Module):
        def __init__(self):
            super(BlackBox, self).__init__()
            self.l = torch.nn.Linear(7, 1, bias=False)
        def forward(self, x):
            return self.l.forward(x)  # Placeholder operation
    #rs_model = BlackBox()


    pipeline = pipelines.make_pipeline(rs_model, None, None)
    #cfg.data.n_points = 100
    #datapoints = data.load_pipeline_data_dict_from_all_sites(cfg.data)
    #datapoint = datapoints[0]
    t = torch.tensor
    datapoint = (
        (
            {
                'sun_H': {'Cc': t([368.0680]), 'IPAR': t([0.]), 'Csl': t([362.6800]), 'ra': t([38.7362]), 'rb': t([19.6884]), 'Ts': t([-5.3135]), 'Pre': t([823.]), 'Ds': t([50.7378]), 'Psi_L': t([0.]), 'Rn': t([-58.5669]), 'QG': t([-68.8654]), 'Vmax': t([25.3091]), 'Psi_sto_50': t([-2.5000]), 'Psi_sto_00': t([-0.5000]), 'CT': t([3.]), 'Ha': t([72.]), 'FI': t([0.0810]), 'Oa': t([210000.]), 'Do': t([800.]), 'a1': t([5.]), 'go': t([0.0100]), 'gmes': t([np.inf]), 'rjv': t([2.1000]), 'DS': t([0.6490]), 'DS_': 0.5},
                'sun_L': {'Cc': t([368.0680]), 'IPAR': t([0.]), 'Csl': t([362.6800]), 'ra': t([38.7362]), 'rb': t([19.6884]), 'Ts': t([-5.3135]), 'Pre': t([823.]), 'Ds': t([50.7378]), 'Psi_L': t([0.]), 'Rn': t([-58.5669]), 'QG': t([-68.8654]), 'Vmax': t([25.3091]), 'Psi_sto_50': t([-2.5000]), 'Psi_sto_00': t([-0.5000]), 'CT': t([3.]), 'Ha': t([72.]), 'FI': t([0.0810]), 'Oa': t([210000.]), 'Do': t([800.]), 'a1': t([5.]), 'go': t([0.0100]), 'gmes': t([np.inf]), 'rjv': t([2.1000]), 'DS': t([0.6490]), 'DS_': 0.5},
                'shd_H': {'Cc': t([368.0680]), 'IPAR': t([0.]), 'Csl': t([362.6800]), 'ra': t([38.7362]), 'rb': t([19.6884]), 'Ts': t([-5.3135]), 'Pre': t([823.]), 'Ds': t([50.7378]), 'Psi_L': t([0.]), 'Rn': t([-58.5669]), 'QG': t([-68.8654]), 'Vmax': t([25.3091]), 'Psi_sto_50': t([-2.5000]), 'Psi_sto_00': t([-0.5000]), 'CT': t([3.]), 'Ha': t([72.]), 'FI': t([0.0810]), 'Oa': t([210000.]), 'Do': t([800.]), 'a1': t([5.]), 'go': t([0.0100]), 'gmes': t([np.inf]), 'rjv': t([2.1000]), 'DS': t([0.6490]), 'DS_': 0.5},
                'shd_L': {'Cc': t([368.0680]), 'IPAR': t([0.]), 'Csl': t([362.6800]), 'ra': t([38.7362]), 'rb': t([19.6884]), 'Ts': t([-5.3135]), 'Pre': t([823.]), 'Ds': t([50.7378]), 'Psi_L': t([0.]), 'Rn': t([-58.5669]), 'QG': t([-68.8654]), 'Vmax': t([25.3091]), 'Psi_sto_50': t([-2.5000]), 'Psi_sto_00': t([-0.5000]), 'CT': t([3.]), 'Ha': t([72.]), 'FI': t([0.0810]), 'Oa': t([210000.]), 'Do': t([800.]), 'a1': t([5.]), 'go': t([0.0100]), 'gmes': t([np.inf]), 'rjv': t([2.1000]), 'DS': t([0.6490]), 'DS_': 0.5},
            }, 
            {'EIn_H': t([0.]), 'EIn_L': t([0.]), 'EG': t([0.0028]), 'ELitter': t([0.]), 'ESN': t([0.]), 'ESN_In': t([0.]), 'EWAT': t([0.]), 'EICE': t([0.]), 'EIn_urb': t([0.]), 'EIn_rock': t([0.])}), 
        t(2.6774)
    )

    print(f"found datapoint {datapoint}")
    x, y = datapoint
    output = pipeline(x)
    graph = torchviz.make_dot(output, params=dict(rs_model.named_parameters()))
    graph.render("torchviz", view=True)

# CONSIDER THIS ONE DONE
def generate_empirical_learning_scatter(config_name):
    # QUEST: Do we simply use random samples as we did usually, or do we try to use actual sample points from sites?
    # The second one seems quite more annoying
    config = OmegaConf.load(f"configs/{config_name}.yaml")
    rs_model = train.train_rs(config)

    rs_sampler = train.rs_sampler()
    with torch.no_grad():
        def point_sampler():
            x = rs_sampler()
            y = train.empirical_model(x)
            return (x, y)
        

        sample_points = train.generate_points(point_sampler, 200)
        path = os.path.join('figures', 'scatter_plots', f'dummy_rs_model_scattering.png')

        plot.fit_plot([rs_model], sample_points, ['our model'], path, f'rs [s/m]',show=True)

# TODO: Do daily aggregates
def generate_Q_LE_timeseries(ordered_data, rs_model):
    # Ordered data is a list of list of points.
    # Inner list of points are from a single site and ordered in time.
    # QUEST: What do we do with missing points ? (NaN, negative Q_LE, etc.)?

    pipeline = pipelines.make_pipeline(rs_model, None, None)
    empirical_pipeline = pipelines.make_pipeline(None, None, None)

    with torch.no_grad():
        for site_points, site_name in ordered_data:
            ys = []
            predicted_ys = []
            empirical_ys = []
            for x,y in site_points:
                ys.append(y)
                predicted_ys.append(pipeline(x))
                empirical_ys.append(empirical_pipeline(x))
            timeseries = [ys, predicted_ys, empirical_ys]
            labels = ['ground truth', 'our pipeline', 'empirical pipeline']
            path = os.path.join('figures', 'timeseries_plots', f'best_model_{site_name}.png')
            plot.time_series_plot(timeseries, labels, path, show=True) 

# TODO: Find interesting x_0 values
def generate_Q_LE_rs_sensitivity_plots():
    rs_range = (10, 3000)

    def Q_LE_Wrapper(x):
        # rs, ra, Rn, QG, Ds, Ts,        
        return differentiable_relations.Q_LE(x[:, 0], x[:, 1], x[:, 2], x[:, 3], x[:, 4], x[:,5])
    
    x_0s = [
        torch.tensor([0., 1., 1., 1., 100., 15.,])
    ]

    x_0_labels = [
        'dummy data'
    ]

    path = os.path.join('figures', 'sensitivity_plots', f'Q_LE_sensitivity_plot.png')
    plot.plot_univariate_slices([Q_LE_Wrapper], x_0s, x_0_labels, 0, 'rs [s/m]', rs_range,
                                1000, ['PM Equation'], 'Q_LE [W/m²]', path=path, show=True)

# TODO: Actually train several models (of course after I'm fixed on a model)
# TODO: Choose interesting datapoints (can take the same as in the normal slices)
# TODO: Also include a table with results for all models (mainly R²?)
def generate_multiple_model_plots(config_names):

    print(f"received config_names={config_names}")

    rs_models = [
        load_model(config_name)
        for config_name in config_names
    ]

    print(f"loaded rs_models={rs_models}")

    for rs_model in rs_models:
        rs_model.eval()

    # scatter : will be hard to see and understand stuff
    # slice : can be interesting
    # timeseries : perhaps a bit less interesting -> because we don't know on what predictors they vary

    x_0s = [
            torch.tensor([-0.2730, 91887.5000, 38.0901, 2.2750, 393.0414, 1000., 14.6980]),
            #torch.tensor([2.510115,   84400,   24.540,  2.2750,  898.362,   800.00,  12.1075])
    ]

    x_0_labels = ['']

    # wrap the rs_model so that it outputs the actual rs
    def wrapper(rs_model, x):
        print(f"in wrapper: got rs_model = {rs_model}")
        x = x * torch.tensor([1e2, 1e-4, 1e-1, 1e1, 1e-1, 1e-2, 1]) # Optional: makes the wrapper take the actual predictor values as inputs
        model_output = rs_model(x)
        rs_small = torch.nn.functional.softplus(model_output)
        return rs_small * 1e2
    
    wrappers = [(lambda x, rs_model=rs_model : wrapper(rs_model, x)) for rs_model in rs_models]
    # wrappers = []
    # for rs_model in rs_models:
        # print(f"computing wrapper for rs_model={rs_model}")
        # l = (lambda x : wrapper(rs_model, x))
        # out = l(x_0s[0].reshape((1, 7)))
        # print(f"immediately: got out={out}")
        # wrappers.append(l)
# 
    # for i,w in enumerate(wrappers):
        # print(f"wrapper {i} : predicting rs={w(x_0s[0].reshape((1, 7)))}")
    
    idx_and_ranges = [
        [1, (82000, 102000), 'Pressure [Pa]'], # Pre plot
        [2, (25, 50), 'Cc [Pa*molCO2/molAIR]'], # Cc plot
        [4, (0, 4000), 'VPD [Pa]'], # Ds plot
        [6, (-10, 40), 'Temperature [°C]'], # Ts plot
    ]

    # for each of the input variables, make a slice plot
    with torch.no_grad():
        for idx, idx_range, idx_label in idx_and_ranges:
            path = os.path.join('figures', 'multiple_models_slices', f'slice_varying_{idx}.png')
            plot.plot_univariate_slices(wrappers, x_0s, x_0_labels, idx, idx_label, idx_range, 100,
                                        config_names, res_label='rs [s/m]',path=path, show=True)



def generate_ctx_decomposition_plots():
    ...
    # NOTE: This one might be tedious to make, and perhaps not the most relevant, 
    # so might drop it if others are enough (or if I run out of time)

# TODO: Replace with the correct data once I choose other sites again
def generate_site_experiment_plots(test_data):
    # scatter plots: interesting to see biases
    # timeseries: perharps harder to interpret
    # slices: hard to interpret. how to choose x0: for sites like the training one, or site like unseen one?
    # now: generate scatter plots, and also give a table with R² scores

    # Load models trained in specific conditions
    hot_rs = load_model('hot_sites_only')
    cold_rs = load_model('cold_sites_only')
    dry_rs = load_model('dry_sites_only')
    wet_rs = load_model('wet_sites_only')

    # Load the model trained with all the sites
    balanced_rs = load_model('best_model')

    hot_rs.eval()
    cold_rs.eval()
    dry_rs.eval()
    wet_rs.eval()
    balanced_rs.eval()

    hot_pipeline = pipelines.make_pipeline(hot_rs, None, None)
    cold_pipeline = pipelines.make_pipeline(cold_rs, None, None)
    dry_pipeline = pipelines.make_pipeline(dry_rs, None, None)
    wet_pipeline = pipelines.make_pipeline(wet_rs, None, None)
    balanced_pipeline = pipelines.make_pipeline(balanced_rs, None, None)

    # Load data from specific sites
    # NOTE: The idea would be to use sites no model was trained on. To avoid suffering from any kind of overfitting,
    # and to focus on the actual climate conditions
    # TODO: That's dumb. better to put maps from site names to data. Because they have intersections, we load sites multiple times
    _, hot_data = load_train_test_data('hot_sites_only')
    _, cold_data = load_train_test_data('cold_sites_only')
    _, dry_data = load_train_test_data('dry_sites_only')
    _, wet_data = load_train_test_data('wet_sites_only')

    hot_data = hot_data[:1000]
    cold_data = cold_data[:1000]
    dry_data = dry_data[:1000]
    wet_data = wet_data[:1000]

    # First, plot the hot, cold and balanced pipeline with cold data
    models_list = [hot_pipeline, cold_pipeline, balanced_pipeline]
    labels = ['hot', 'cold', 'full']
    path = os.path.join('figures', 'site_experiment_scatters', f'cold_on_hot_scatter.png')
    with torch.no_grad():
        plot.fit_plot(models_list, cold_data, labels, path, plot_range=(0, 600), unit=f'Q_LE [W/m²]',show=True)

    # Then, plot the hot, cold and balanced pipeline with hot data
    path = os.path.join('figures', 'site_experiment_scatters', f'hot_on_cold_scatter.png')
    with torch.no_grad():
        plot.fit_plot(models_list, hot_data, labels, path, plot_range=(0, 600), unit=f'Q_LE [W/m²]',show=True)

    # Now: plot the dry, wet and balanced pipeline with wet data
    models_list = [dry_pipeline, wet_pipeline, balanced_pipeline]
    labels = ['dry', 'wet', 'full']
    path = os.path.join('figures', 'site_experiment_scatters', f'wet_on_dry_scatter.png')
    with torch.no_grad():
        plot.fit_plot(models_list, wet_data, labels, path, plot_range=(0, 600), unit=f'Q_LE [W/m²]',show=True)

    # Then, plot the dry, wet and balanced pipeline with dry data
    path = os.path.join('figures', 'site_experiment_scatters', f'dry_on_wet_scatter.png')
    with torch.no_grad():
        plot.fit_plot(models_list, dry_data, labels, path, plot_range=(0, 600), unit=f'Q_LE [W/m²]',show=True)


# Objects and data needed for plotting
#rs_model = load_model('best_model')
#train_data, test_data = load_train_test_data('best_model')
#train_data, test_data = load_train_test_data('default')
#ordered_data = load_ordered_data('best_model')
multiple_model_config_names = [
    'best_model',
    'cold_sites_only',
    'hot_sites_only',
    'dry_sites_only',
    'default'
]

# Call the desired functions generating plots
#generate_hp_tuning_loss_plot()
#generate_torch_viz_graph()
#generate_rs_slice_plots(rs_model)
#generate_scatter_plots(test_data, rs_model)
#generate_Q_LE_timeseries(ordered_data, rs_model)
#generate_empirical_learning_scatter('best_model')
#generate_Q_LE_rs_sensitivity_plots()
#generate_multiple_model_plots(multiple_model_config_names)
#generate_Q_LE_slice_plots(rs_model)
generate_site_experiment_plots(None)


