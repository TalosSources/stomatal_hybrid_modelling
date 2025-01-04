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

import torchviz

import plot
import pipelines
import models
import data
import differentiable_relations

import train

def load_model(config_name): # TODO: yaml and config names might be different, fix that (hot_sites_only, hot_sites_only_experiment.yaml)
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

    coefficients = [ # TODO: Replace with actual coefficients
        0.25,
        0.36,
        0.78
    ]

    # TODO: Latex table with coefficients? or perhaps another tool for the table

    plot.plot_losses(losses, labels, coefficients=coefficients, minmax=False, smoothing=0.01)

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

        # Define what to plot TODO: Check units
        # TODO: Make a common plot with all chosen variables? plt.subplot... requires factoring everything...
        idx_and_ranges = [
            [1, (82000, 102000), 'Pressure [Pa]'], # Pre plot
            [2, (25, 50), 'Cc'], # Cc plot
            [4, (0, 2500), 'Ds [Pa]'], # Ds plot
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

def generate_Q_LE_slice_plots():
    # TODO
    # Idea: work with lambda. have a lambda closure that contains 'hardcoded' (using fetched data) predictors for all values except
    # the ones that we vary
    # load specific(s) datapoint(s)
    # I'll vary Ds, Ts, Cc, and Pre.
    # An=[-0.2730], Pre=[91887.5000], Cc=[38.0901], GAM=[2.2750], Ds=[393.0414], Do=[1000.], Ts=[14.6980]
    with torch.no_grad():
        #x_0 = torch.tensor([-0.2730, 91887.5000, 38.0901, 2.2750, 393.0414, 1000., 14.6980])
        x_0 = torch.tensor([2.510115,   84400,   24.540,  2.2750,  898.362,   800.00,  12.1075])

        # load the model weights
        results_path = os.path.join('results', 'default')
        rs_model_state_dict = torch.load(os.path.join(results_path, "model_weights.pt"))
        cfg = OmegaConf.load("configs/default.yaml")
        rs_model = models.rs_model_from_config(cfg.model)
        rs_model.load_state_dict(rs_model_state_dict)

        # wrap the rs_model so that it outputs the actual rs
        def wrapper(x):
            x = x * torch.tensor([1e2, 1e-4, 1e-1, 1e1, 1e-1, 1e-2, 1]) # Optional: makes the wrapper take the actual predictor values as inputs
            model_output = rs_model(x)
            rs_small = torch.nn.functional.softplus(model_output)
            return rs_small * 1e2
        
        pipeline = pipelines.make_pipeline(rs_model)
        
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

        # Define what to plot TODO: Check units
        # TODO: Make a common plot with all chosen variables? plt.subplot... requires factoring everything...
        idx_and_ranges = [
            [1, (82000, 102000), 'Pressure [Pa]'], # Pre plot
            [2, (25, 50), 'Cc'], # Cc plot
            [4, (0, 2500), 'Ds [Pa]'], # Ds plot
            [6, (-10, 40), 'Temperature [°C]'], # Ts plot
        ]

        # we might add more models on the same plot
        labels = ['rs_model (ours)', 'empirical model']
        models_list = [wrapper, empirical_model]

        # for each of the input variables, make a slice plot
        for idx, idx_range, idx_label in idx_and_ranges:
            path = os.path.join('figures', 'rs_slices', f'slice_varying_{idx}.png')
            plot.plot_univariate_slices(models_list, x_0, idx, idx_label, idx_range, 100, labels, path=path, show=True)

def generate_scatter_plots(test_data, rs_model):
    # TODO: Data, rs_model, that kind of things could be factorized, to avoid doing it that many times
    # TODO: Do daily aggregates. Idea: log in global_preds, the i index of the timestep. 
    # then, write a function that, given data from some site (in sequential order? perhaps not needed),
    # groups the data in  
    
    test_data = test_data[:5000] # take only a subset of the data. for performances. TODO: Use ctx batches in the plot method (or here?), for better perf.

    # make the pipeline
    pipeline = pipelines.make_pipeline(rs_model, None, None)
    print(f"finished making pipeline")

    empirical_pipeline = pipelines.make_pipeline(None, None, None)

    models_list = [pipeline, empirical_pipeline]
    labels = ['our pipeline', 'empirical pipeline']

    path = os.path.join('figures', 'scatter_plots', f'simple_best_model_scattering.png')

    with torch.no_grad():
        plot.fit_plot(models_list, test_data, labels, path, plot_range=(0, 600), show=True)

def generate_torch_viz_graph():
    # TODO: For now it only shows the model architecture, and some long line.
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
    pipeline = pipelines.make_pipeline(rs_model, None, None)
    datapoints = data.load_pipeline_data_dict_from_all_sites(cfg.data)
    datapoint = datapoints[0]
    x, y = datapoint
    output = pipeline(x)
    graph = torchviz.make_dot(output, params=dict(rs_model.named_parameters()))
    graph.render("torchviz", view=True)

def generate_empirical_learning_scatter(config_name):
    # CONSIDER THIS ONE DONE
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

def generate_Q_LE_timeseries(ordered_data, rs_model):
    # Ordered data is a list of list of points.
    # Inner list of points are from a single site and ordered in time.
    # TODO / QUEST: What do we do with missing points ? (NaN, negative Q_LE, etc.)?

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

def generate_Q_LE_rs_sensitivity_plots():
    rs_range = (10, 3000)

    def Q_LE_Wrapper(x):
        # rs, ra, Rn, QG, Ds, Ts,        
        return differentiable_relations.Q_LE(x[:, 0], x[:, 1], x[:, 2], x[:, 3], x[:, 4], x[:,5])
    
    x_0s = [
        torch.tensor([0., 1., 1., 1., 100., 15.,]) # TODO: Replace with actual values (and several datapoints from several datasets?)
    ]

    x_0_labels = [
        'dummy data'
    ]

    path = os.path.join('figures', 'sensitivity_plots', f'Q_LE_sensitivity_plot.png')
    plot.plot_univariate_slices([Q_LE_Wrapper], x_0s, x_0_labels, 0, 'rs [s/m]', rs_range,
                                1000, ['PM Equation'], 'Q_LE [W/m²]', path=path, show=True)

def generate_multiple_model_plots():
    ... # TODO: Choose which plot to make (timeserie, scatter, slice)

def generate_ctx_decomposition_plots():
    ... # TODO
    # NOTE: This one might be tedious to make, and perhaps not the most relevant, 
    # so might drop it if others are enough (or if I run out of time)

def generate_site_experiment_plots():
    ...


# Objects and data needed for plotting
#rs_model = load_model('best_model')
#train_data, test_data = load_train_test_data('best_model')
#ordered_data = load_ordered_data('best_model')

# Call the desired functions generating plots
#generate_hp_tuning_loss_plot()
#generate_torch_viz_graph()
#generate_rs_slice_plots(rs_model)
#generate_scatter_plots(test_data, rs_model)
#generate_Q_LE_timeseries(ordered_data, rs_model)
#generate_empirical_learning_scatter('best_model')
generate_Q_LE_rs_sensitivity_plots()


