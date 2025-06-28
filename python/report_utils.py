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
import eval

import train


def load_model(config_name):
    # load the model weights
    results_path = os.path.join("results", config_name)
    rs_model_state_dict = torch.load(os.path.join(results_path, "model_weights.pt"))
    cfg = OmegaConf.load(f"configs/{config_name}.yaml")
    rs_model = models.rs_model_from_config(cfg.model)
    rs_model.load_state_dict(rs_model_state_dict)
    return rs_model


def load_train_test_data(config_name):
    config = OmegaConf.load(f"configs/{config_name}.yaml")
    all_data = data.load_pipeline_data_dict_from_all_sites(
        config.data, verbose=True, output_keys=["LE"]
    )  # NOTE: LE_CORR?
    train_data, test_data = train_test_split(
        all_data,
        test_size=config.train.test_split,
        random_state=config.train.seed,
        shuffle=True,
    )
    return train_data, test_data


def load_ordered_data(config_name):
    config = OmegaConf.load(f"configs/{config_name}.yaml")
    site_names = config.data.sites
    base_path = os.path.expanduser(config.data.base_path)
    load_ordered_data_from_sites(site_names, base_path)


def load_ordered_data_from_sites(site_names, base_path):
    combined_site_data = []
    for site_name in site_names:
        pred_path = os.path.join(base_path, f"Results_{site_name}.mat")
        obs_path = os.path.join(base_path, f"Res_{site_name}.mat")
        exclude_negative_outputs = True
        site_data = []
        n_points = 25000  # 1000 corresponds to about 3 weeks NOTE: Where should we choose this?
        data.load_pipeline_data_dict_single_site(
            pred_path,
            obs_path,
            site_data,
            nPoints=n_points,
            exclude_negative_outputs=exclude_negative_outputs,
            verbose=True,
        )
        if len(site_data) > 0:
            combined_site_data.append((site_data, site_name))
    return combined_site_data


def generate_hp_tuning_loss_plot():
    results_path = "results/lr_wd_tuning"
    filenames = [
        "losses_lr-0.001_weight_decay-1e-06_n_hidden-4_hidden_size-128_batch_norm-False_.npy",
        "losses_lr-0.0001_weight_decay-1e-06_n_hidden-4_hidden_size-128_batch_norm-False_.npy",
        "losses_lr-1e-05_weight_decay-0_n_hidden-4_hidden_size-128_batch_norm-False_.npy",
    ]

    losses = [np.load(os.path.join(results_path, filename)) for filename in filenames]

    labels = [
        "lr=1e-3, λ=1e-6",
        "lr=1e-4, λ=1e-6",
        "lr=1e-5, λ=0",
    ]

    coefficients = [
        0.638,
        0.622,
        0.537,
    ]

    plot.plot_losses(
        losses, labels, coefficients=coefficients, minmax=False, smoothing=0.03
    )


# CONSIDER THIS DONE
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
        # plots: one in cold, one in hot, one in temperate arid, one in temperate wet
        # cold/hot: use temperatures, and take cold/hot sites
        # arid/wet: use vpd. low vpd means transpiration could happen: some notion of "dry"
        # cold+humid: An=tensor([1.2889]), Pre=tensor([83600.]), Cc=tensor([26.8077]), GAM=tensor([0.8524]), Ds=tensor([185.5382]), Do=tensor([800.]), Ts=tensor([-0.5985]) (CH-Dav)
        # hot+dry: An=tensor([-1.3149]), Pre=tensor([97125.]), Cc=tensor([51.0093]), GAM=tensor([7.0099]), Ds=tensor([4388.1724]), Do=tensor([1500.]), Ts=tensor([35.5950]) (AU-Stp)
        # cold+dry An=tensor([-0.2617]), Pre=tensor([101659.5000]), Cc=tensor([42.6090]), GAM=tensor([1.4971]), Ds=tensor([395.3699]), Do=tensor([1000.]), Ts=tensor([5.5390])(ES-Amo)
        # hot+humid An=tensor([-0.3253]), Pre=tensor([100600.]), Cc=tensor([41.7964]), GAM=tensor([4.0832]), Ds=tensor([222.2412]), Do=tensor([1000.]), Ts=tensor([23.9900]) (GF-Guy)
        x_0s = [
            torch.tensor([1.2889, 83600.0, 26.8077, 0.8524, 185.5382, 800.0, -0.5985]),
            torch.tensor(
                [-1.3149, 97125.0, 51.0093, 7.0099, 4388.1724, 1500.0, 35.5950]
            ),
            torch.tensor(
                [-0.2617, 101659.5000, 42.6090, 1.4971, 395.3699, 1000.0, 5.5390]
            ),
            torch.tensor(
                [-0.3253, 100600.0, 41.7964, 4.0832, 222.2412, 1000.0, 23.9900]
            ),
        ]

        x_0_labels = [
            "cold, humid (CH-Dav)",
            "hot, dry (AU-Stp)",
            "cold, dry (ES-Amo)",
            "hot, humid (GF-Guy)",
        ]

        # wrap the rs_model so that it outputs the actual rs
        def wrapper(x):
            x = x * torch.tensor(
                [1e2, 1e-4, 1e-1, 1e1, 1e-1, 1e-2, 1]
            )  # Optional: makes the wrapper take the actual predictor values as inputs
            model_output = rs_model(x)
            rs_small = torch.nn.functional.softplus(model_output)
            return rs_small * 1e2

        def empirical_model(x):
            a1 = 5
            go = 10000
            An = x[:, 0]
            Pre = x[:, 1]
            Cc = x[:, 2]
            GAM = x[:, 3]
            Ds = x[:, 4]
            Do = x[:, 5]
            Ts = x[:, 6]
            gsCO2 = go + a1 * An * Pre / (
                (Cc - GAM) * (1 + Ds / Do)
            )  ###  [umolCO2 / s m^2] -- Stomatal Conductance
            gsCO2 = torch.max(gsCO2, torch.tensor(go))
            rsCO2 = 1 / gsCO2  ### [ s m^2 / umolCO2 ] Stomatal resistence or Canopy
            rsH20 = (rsCO2 / 1.64) * (
                1e6
            )  ### [ s m^2 / molH20 ] Stomatal resistence or canopy
            rs = (
                rsH20 * (273.15 * Pre) / (0.0224 * (Ts + 273.15) * 101325)
            )  ## [s/m]  Stomatal resistence or Canopy [convert stomatal resistence in terms of water volumes into another unit]
            return rs

        # Define what to plot
        idx_and_ranges = [
            [1, (82000, 102000), "Pressure [Pa]"],  # Pre plot
            [2, (25, 50), "Cc [Pa*molCO2/molAIR]"],  # Cc plot
            [4, (0, 4000), "VPD [Pa]"],  # Ds plot
            [6, (-10, 40), "Temperature [°C]"],  # Ts plot
        ]

        # we might add more models on the same plot
        labels = ["Hybrid model (ours)", "Empirical model"]
        models_list = [wrapper, empirical_model]

        path = os.path.join("figures", "rs_slices", f"combined_slices.png")
        plot.plot_univariate_slices_subplots(
            models_list,
            x_0s,
            x_0_labels,
            idx_and_ranges,
            100,
            labels,
            res_label="rs [s/m]",
            path=path,
            show=True,
        )


# CONSIDER THIS DONE
def generate_Q_LE_slice_plots(rs_model):
    with torch.no_grad():
        t = torch.tensor

        # Hot, Humid (GF-Guy)
        hot_wet_point = (
            (
                {
                    "sun_H": {
                        "Cc": t([289.3110]),
                        "IPAR": t([90.7115]),
                        "Csl": t([376.3200]),
                        "ra": t([8.8555]),
                        "rb": t([51.0891]),
                        "Ts": t([28.1250]),
                        "Pre": t([1007.0]),
                        "Ds": t([967.4290]),
                        "Psi_L": t([-0.0663]),
                        "Rn": t([475.5663]),
                        "QG": t([127.6617]),
                        "Vmax": t([37.7132]),
                        "Psi_sto_50": t([-1.7000]),
                        "Psi_sto_00": t([-0.9000]),
                        "CT": t([3.0]),
                        "Ha": t([72.0]),
                        "FI": t([0.0810]),
                        "Oa": t([210000.0]),
                        "Do": t([1000.0]),
                        "a1": t([6.0]),
                        "go": t([0.0100]),
                        "gmes": t([np.inf]),
                        "rjv": t([2.2000]),
                        "DS": t([0.6490]),
                    },
                    "shd_H": {
                        "Cc": t([298.4520]),
                        "IPAR": t([35.0558]),
                        "Csl": t([376.3200]),
                        "ra": t([8.8555]),
                        "rb": t([51.0891]),
                        "Ts": t([28.1250]),
                        "Pre": t([1007.0]),
                        "Ds": t([967.4290]),
                        "Psi_L": t([-0.0663]),
                        "Rn": t([475.5663]),
                        "QG": t([127.6617]),
                        "Vmax": t([23.1776]),
                        "Psi_sto_50": t([-1.7000]),
                        "Psi_sto_00": t([-0.9000]),
                        "CT": t([3.0]),
                        "Ha": t([72.0]),
                        "FI": t([0.0810]),
                        "Oa": t([210000.0]),
                        "Do": t([1000.0]),
                        "a1": t([6.0]),
                        "go": t([0.0100]),
                        "gmes": t([np.inf]),
                        "rjv": t([2.2000]),
                        "DS": t([0.6490]),
                    },
                },
                {
                    "EIn_H": t([0.1231]),
                    "EIn_L": t([0.0]),
                    "EG": t([0.0432]),
                    "ELitter": t([0.0]),
                    "ESN": t([0.0]),
                    "ESN_In": t([0.0]),
                    "EWAT": t([0.0]),
                    "EICE": t([0.0]),
                    "EIn_urb": t([0.0]),
                    "EIn_rock": t([0.0]),
                },
            ),
            t(338.4500),
        )

        # Cold, Humid (CH-Dav)
        cold_wet_point = (
            (
                {
                    "sun_H": {
                        "Cc": t([319.2088]),
                        "IPAR": t([40.6829]),
                        "Csl": t([362.6200]),
                        "ra": t([61.0868]),
                        "rb": t([24.1148]),
                        "Ts": t([-4.3020]),
                        "Pre": t([829.9150]),
                        "Ds": t([92.6856]),
                        "Psi_L": t([-0.0258]),
                        "Rn": t([9.2711]),
                        "QG": t([-0.6741]),
                        "Vmax": t([35.1236]),
                        "Psi_sto_50": t([-2.5000]),
                        "Psi_sto_00": t([-0.5000]),
                        "CT": t([3.0]),
                        "Ha": t([72.0]),
                        "FI": t([0.0810]),
                        "Oa": t([210000.0]),
                        "Do": t([800.0]),
                        "a1": t([5.0]),
                        "go": t([0.0100]),
                        "gmes": t([np.inf]),
                        "rjv": t([2.1000]),
                        "DS": t([0.6490]),
                    },
                    "shd_H": {
                        "Cc": t([326.8221]),
                        "IPAR": t([16.2409]),
                        "Csl": t([362.6200]),
                        "ra": t([61.0868]),
                        "rb": t([24.1148]),
                        "Ts": t([-4.3020]),
                        "Pre": t([829.9150]),
                        "Ds": t([92.6856]),
                        "Psi_L": t([-0.0258]),
                        "Rn": t([9.2711]),
                        "QG": t([-0.6741]),
                        "Vmax": t([23.2060]),
                        "Psi_sto_50": t([-2.5000]),
                        "Psi_sto_00": t([-0.5000]),
                        "CT": t([3.0]),
                        "Ha": t([72.0]),
                        "FI": t([0.0810]),
                        "Oa": t([210000.0]),
                        "Do": t([800.0]),
                        "a1": t([5.0]),
                        "go": t([0.0100]),
                        "gmes": t([np.inf]),
                        "rjv": t([2.1000]),
                        "DS": t([0.6490]),
                    },
                },
                {
                    "EIn_H": t([0.0]),
                    "EIn_L": t([0.0]),
                    "EG": t([0.0]),
                    "ELitter": t([0.0]),
                    "ESN": t([0.0083]),
                    "ESN_In": t([0.0025]),
                    "EWAT": t([0.0]),
                    "EICE": t([0.0]),
                    "EIn_urb": t([0.0]),
                    "EIn_rock": t([0.0]),
                },
            ),
            t(25.3235),
        )

        # Hot, Dry (AU-STP)
        hot_dry_point = (
            (
                {
                    "sun_L": {
                        "Cc": t([223.7393]),
                        "IPAR": t([247.4754]),
                        "Csl": t([386.0600]),
                        "ra": t([107.2802]),
                        "rb": t([23.3581]),
                        "Ts": t([30.8950]),
                        "Pre": t([974.6000]),
                        "Ds": t([1786.9075]),
                        "Psi_L": t([-0.3937]),
                        "Rn": t([596.9096]),
                        "QG": t([157.3113]),
                        "Vmax": t([30.4333]),
                        "Psi_sto_50": t([-2.0]),
                        "Psi_sto_00": t([-0.3000]),
                        "CT": t([4.0]),
                        "Ha": t([72.0]),
                        "FI": t([0.0400]),
                        "Oa": t([210000.0]),
                        "Do": t([1500.0]),
                        "a1": t([6.0]),
                        "go": t([0.0100]),
                        "gmes": t([np.inf]),
                        "rjv": t([1.9000]),
                        "DS": t([0.6490]),
                    },
                    "shd_L": {
                        "Cc": t([237.9602]),
                        "IPAR": t([174.7201]),
                        "Csl": t([386.0600]),
                        "ra": t([107.2802]),
                        "rb": t([23.3581]),
                        "Ts": t([30.8950]),
                        "Pre": t([974.6000]),
                        "Ds": t([1786.9075]),
                        "Psi_L": t([-0.3937]),
                        "Rn": t([596.9096]),
                        "QG": t([157.3113]),
                        "Vmax": t([29.8925]),
                        "Psi_sto_50": t([-2.0]),
                        "Psi_sto_00": t([-0.3000]),
                        "CT": t([4.0]),
                        "Ha": t([72.0]),
                        "FI": t([0.0400]),
                        "Oa": t([210000.0]),
                        "Do": t([1500.0]),
                        "a1": t([6.0]),
                        "go": t([0.0100]),
                        "gmes": t([np.inf]),
                        "rjv": t([1.9000]),
                        "DS": t([0.6490]),
                    },
                },
                {
                    "EIn_H": t([0.0]),
                    "EIn_L": t([0.0]),
                    "EG": t([0.0631]),
                    "ELitter": t([0.0]),
                    "ESN": t([0.0]),
                    "ESN_In": t([0.0]),
                    "EWAT": t([0.0]),
                    "EICE": t([0.0]),
                    "EIn_urb": t([0.0]),
                    "EIn_rock": t([0.0]),
                },
            ),
            t(293.0),
        )

        x_0s = [
            (torch.tensor([1007.0, 967.4290, 28.1250]), hot_wet_point),
            (torch.tensor([829.9150, 92.6856, -4.3020]), cold_wet_point),
            (torch.tensor([974.6000, 1786.9075, 30.8950]), hot_dry_point),
        ]
        x_0_labels = ["hot, humid", "cold, humid", "hot, dry"]

        pipeline = pipelines.make_pipeline(rs_model, None, None)
        empirical_pipeline = pipelines.make_pipeline(None, None, None)

        steps = 100

        def wrapper(pip, x, datapoint):
            # x is [Pre, Ds, Ts]
            datapoint_2 = deepcopy(datapoint)
            p, _ = datapoint_2
            ctx_preds, _ = p
            for _, pred in ctx_preds.items():  # NOTE: Does it vary datapoint_2?
                pred["Pre"] = x[:, 0]
                pred["Ds"] = x[:, 1]
                pred["Ts"] = x[:, 2]
                for k, v in pred.items():
                    if k not in ["Pre", "Ds", "Ts"]:
                        pred[k] = pred[k].expand(steps)
            return pip(p)

        our_wrapper = lambda x, datapoint: wrapper(pipeline, x, datapoint)
        empirical_wrapper = lambda x, datapoint: wrapper(
            empirical_pipeline, x, datapoint
        )

        idx_and_ranges = [
            [0, (820, 1020), "Pressure [Pa]"],  # Pre plot
            [1, (0, 4000), "VPD [Pa]"],  # Ds plot
            [2, (-10, 40), "Temperature [°C]"],  # Ts plot
        ]

        # we might add more models on the same plot
        labels = ["Hybrid Model (ours)", "Empirical Model"]
        models_list = [our_wrapper, empirical_wrapper]

        path = os.path.join("figures", "Q_LE_slices", f"combined_slices.png")
        plot.plot_univariate_slices_subplots(
            models_list,
            x_0s,
            x_0_labels,
            idx_and_ranges,
            steps,
            labels,
            res_label="Q_LE [W/m²]",
            path=path,
            x_0_suppl=True,
            show=True,
        )


# NOTE: We could do daily aggregates. Idea: log in global_preds, the i index of the timestep.
def generate_scatter_plots(test_data, rs_model):
    # then, write a function that, given data from some site (in sequential order? perhaps not needed),
    # groups the data in

    test_data = test_data[
        :5000
    ]  # take only a subset of the data. for performances. OPT: Use ctx batches in the plot method (or here?), for better perf because it's so slow now.

    # make the pipeline
    pipeline = pipelines.make_pipeline(rs_model, None, None)

    empirical_pipeline = pipelines.make_pipeline(None, None, None)

    models_list = [pipeline, empirical_pipeline]
    labels = ["Hybrid model", "Empirical model"]

    path = os.path.join("figures", "scatter_plots", f"simple_best_model_scattering.png")

    with torch.no_grad():
        plot.fit_plot(
            models_list,
            test_data,
            labels,
            path,
            plot_range=(0, 600),
            unit=f"Q_LE [W/m²]",
            show=True,
        )


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

    # rs_model = BlackBox()

    pipeline = pipelines.make_pipeline(rs_model, None, None)
    # cfg.data.n_points = 100
    # datapoints = data.load_pipeline_data_dict_from_all_sites(cfg.data)
    # datapoint = datapoints[0]
    t = torch.tensor
    datapoint = (
        (
            {
                "sun_H": {
                    "Cc": t([368.0680]),
                    "IPAR": t([0.0]),
                    "Csl": t([362.6800]),
                    "ra": t([38.7362]),
                    "rb": t([19.6884]),
                    "Ts": t([-5.3135]),
                    "Pre": t([823.0]),
                    "Ds": t([50.7378]),
                    "Psi_L": t([0.0]),
                    "Rn": t([-58.5669]),
                    "QG": t([-68.8654]),
                    "Vmax": t([25.3091]),
                    "Psi_sto_50": t([-2.5000]),
                    "Psi_sto_00": t([-0.5000]),
                    "CT": t([3.0]),
                    "Ha": t([72.0]),
                    "FI": t([0.0810]),
                    "Oa": t([210000.0]),
                    "Do": t([800.0]),
                    "a1": t([5.0]),
                    "go": t([0.0100]),
                    "gmes": t([np.inf]),
                    "rjv": t([2.1000]),
                    "DS": t([0.6490]),
                    "DS_": 0.5,
                },
                "sun_L": {
                    "Cc": t([368.0680]),
                    "IPAR": t([0.0]),
                    "Csl": t([362.6800]),
                    "ra": t([38.7362]),
                    "rb": t([19.6884]),
                    "Ts": t([-5.3135]),
                    "Pre": t([823.0]),
                    "Ds": t([50.7378]),
                    "Psi_L": t([0.0]),
                    "Rn": t([-58.5669]),
                    "QG": t([-68.8654]),
                    "Vmax": t([25.3091]),
                    "Psi_sto_50": t([-2.5000]),
                    "Psi_sto_00": t([-0.5000]),
                    "CT": t([3.0]),
                    "Ha": t([72.0]),
                    "FI": t([0.0810]),
                    "Oa": t([210000.0]),
                    "Do": t([800.0]),
                    "a1": t([5.0]),
                    "go": t([0.0100]),
                    "gmes": t([np.inf]),
                    "rjv": t([2.1000]),
                    "DS": t([0.6490]),
                    "DS_": 0.5,
                },
                "shd_H": {
                    "Cc": t([368.0680]),
                    "IPAR": t([0.0]),
                    "Csl": t([362.6800]),
                    "ra": t([38.7362]),
                    "rb": t([19.6884]),
                    "Ts": t([-5.3135]),
                    "Pre": t([823.0]),
                    "Ds": t([50.7378]),
                    "Psi_L": t([0.0]),
                    "Rn": t([-58.5669]),
                    "QG": t([-68.8654]),
                    "Vmax": t([25.3091]),
                    "Psi_sto_50": t([-2.5000]),
                    "Psi_sto_00": t([-0.5000]),
                    "CT": t([3.0]),
                    "Ha": t([72.0]),
                    "FI": t([0.0810]),
                    "Oa": t([210000.0]),
                    "Do": t([800.0]),
                    "a1": t([5.0]),
                    "go": t([0.0100]),
                    "gmes": t([np.inf]),
                    "rjv": t([2.1000]),
                    "DS": t([0.6490]),
                    "DS_": 0.5,
                },
                "shd_L": {
                    "Cc": t([368.0680]),
                    "IPAR": t([0.0]),
                    "Csl": t([362.6800]),
                    "ra": t([38.7362]),
                    "rb": t([19.6884]),
                    "Ts": t([-5.3135]),
                    "Pre": t([823.0]),
                    "Ds": t([50.7378]),
                    "Psi_L": t([0.0]),
                    "Rn": t([-58.5669]),
                    "QG": t([-68.8654]),
                    "Vmax": t([25.3091]),
                    "Psi_sto_50": t([-2.5000]),
                    "Psi_sto_00": t([-0.5000]),
                    "CT": t([3.0]),
                    "Ha": t([72.0]),
                    "FI": t([0.0810]),
                    "Oa": t([210000.0]),
                    "Do": t([800.0]),
                    "a1": t([5.0]),
                    "go": t([0.0100]),
                    "gmes": t([np.inf]),
                    "rjv": t([2.1000]),
                    "DS": t([0.6490]),
                    "DS_": 0.5,
                },
            },
            {
                "EIn_H": t([0.0]),
                "EIn_L": t([0.0]),
                "EG": t([0.0028]),
                "ELitter": t([0.0]),
                "ESN": t([0.0]),
                "ESN_In": t([0.0]),
                "EWAT": t([0.0]),
                "EICE": t([0.0]),
                "EIn_urb": t([0.0]),
                "EIn_rock": t([0.0]),
            },
        ),
        t(2.6774),
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

        sample_points = train.generate_points(point_sampler, 1000)
        coeff = eval.simple_coefficient_determination(rs_model, sample_points)
        print(f"found coeff={coeff}")
        path = os.path.join(
            "figures", "scatter_plots", f"dummy_rs_model_scattering.svg"
        )

        plot.fit_plot(
            [rs_model], sample_points, ["our model"], path, f"rs [s/m]", show=True
        )


# CONSIDER THIS DONE
def generate_Q_LE_timeseries(ordered_data, ranges, rs_model):
    # Ordered data is a list of list of points.
    # Inner list of points are from a single site and ordered in time.
    # QUEST: What do we do with missing points ? (NaN, negative Q_LE, etc.)?

    # interesting intervals by site:
    # CN-Du2: 1500-1800
    # US-Me2: 2400-2700?
    # CH-Dav: very variable. maybe smtg like 6400-6700 to show different intensities
    # ES-Amo: Anywhere is good, the empirical model so drastically overestimates
    # Au-Wac: bad throughout, both the empirical and us underestimate, but we do a little more
    # Au-Stp: Anywhere is good, same as ES-Amo but even more
    # Gf-Guy: Most places, empirical underestimates, we do better.

    pipeline = pipelines.make_pipeline(rs_model, None, None)
    empirical_pipeline = pipelines.make_pipeline(None, None, None)

    linestyles = ["solid", "solid", "solid"]
    widths = [1.5, 1.8, 1.5]
    alphas = [0.6, 0.8, 0.6]

    with torch.no_grad():
        for (site_points, site_name), range in zip(ordered_data, ranges):
            ys = []
            predicted_ys = []
            empirical_ys = []
            for x, y in site_points[range[0] : range[1]]:
                ys.append(y)
                predicted_ys.append(pipeline(x))
                empirical_ys.append(empirical_pipeline(x))
                # predicted_ys.append(y)
                # empirical_ys.append(y)
            timeseries = [ys, predicted_ys, empirical_ys]
            labels = ["Observed", "Hybrid Model", "Empirical Model"]
            path = os.path.join(
                "figures", "timeseries_plots", f"best_model_{site_name}.svg"
            )
            plot.time_series_plot(
                timeseries,
                labels,
                path,
                linestyles=linestyles,
                widths=widths,
                alphas=alphas,
                x_label="Timestep",
                y_label="Q_LE [W/m²]",
                show=True,
            )


def generate_multiple_model_plots(rs_models):
    for rs_model in rs_models:
        rs_model.eval()

    # scatter : will be hard to see and understand stuff
    # slice : can be interesting
    # timeseries : perhaps a bit less interesting -> because we don't know on what predictors they vary

    x_0s = [
        torch.tensor([1.2889, 83600.0, 26.8077, 0.8524, 185.5382, 800.0, -0.5985]),
        torch.tensor([-1.3149, 97125.0, 51.0093, 7.0099, 4388.1724, 1500.0, 35.5950]),
        torch.tensor([-0.2617, 101659.5000, 42.6090, 1.4971, 395.3699, 1000.0, 5.5390]),
        torch.tensor([-0.3253, 100600.0, 41.7964, 4.0832, 222.2412, 1000.0, 23.9900]),
    ]

    x_0_labels = [
        "cold, humid (CH-Dav)",
        "hot, dry (AU-Stp)",
        "cold, dry (ES-Amo)",
        "hot, humid (GF-Guy)",
    ]

    # wrap the rs_model so that it outputs the actual rs
    def wrapper(rs_model, x):
        x = x * torch.tensor(
            [1e2, 1e-4, 1e-1, 1e1, 1e-1, 1e-2, 1]
        )  # Optional: makes the wrapper take the actual predictor values as inputs
        model_output = rs_model(x)
        rs_small = torch.nn.functional.softplus(model_output)
        return rs_small * 1e2

    wrappers = [
        (lambda x, rs_model=rs_model: wrapper(rs_model, x)) for rs_model in rs_models
    ]
    labels = [f"Model {i + 1}" for i in range(len(rs_models))]

    idx_and_ranges = [
        [1, (82000, 102000), "Pressure [Pa]"],  # Pre plot
        [2, (25, 50), "Cc [Pa*molCO2/molAIR]"],  # Cc plot
        [4, (0, 4000), "VPD [Pa]"],  # Ds plot
        [6, (-10, 40), "Temperature [°C]"],  # Ts plot
    ]

    # for each of the input variables, make a slice plot
    with torch.no_grad():
        path = os.path.join(
            "figures", "multiple_models_slices", f"combined_multimodel_slices.svg"
        )
        plot.plot_univariate_slices_subplots(
            wrappers,
            x_0s,
            x_0_labels,
            idx_and_ranges,
            100,
            labels,
            res_label="rs [s/m]",
            path=path,
            show=True,
        )


# CONSIDER THIS DONE
def generate_site_experiment_plots(
    hot_rs, cold_rs, dry_rs, wet_rs, balanced_rs, test_ordered_data, seed=24
):
    # scatter plots: interesting to see biases
    # timeseries: perharps harder to interpret
    # slices: hard to interpret. how to choose x0: for sites like the training one, or site like unseen one?
    # now: generate scatter plots, and also give a table with R² scores

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
    dry_data, _ = test_ordered_data[0]
    wet_data, _ = test_ordered_data[1]
    hot_data, _ = test_ordered_data[2]
    cold_data, _ = test_ordered_data[3]

    trunc = 2000
    gen = np.random.default_rng(seed)

    dry_data = np.array(dry_data, dtype=object)[
        gen.choice(
            len(dry_data),
            trunc,
            replace=False,
        )
    ]
    wet_data = np.array(wet_data, dtype=object)[
        gen.choice(
            len(wet_data),
            trunc,
            replace=False,
        )
    ]
    hot_data = np.array(hot_data, dtype=object)[
        gen.choice(
            len(hot_data),
            trunc,
            replace=False,
        )
    ]
    cold_data = np.array(cold_data, dtype=object)[
        gen.choice(
            len(cold_data),
            trunc,
            replace=False,
        )
    ]

    # First, plot the hot, cold and balanced pipeline with cold data
    models_list = [hot_pipeline, cold_pipeline, balanced_pipeline]
    labels = ["Hot Model", "Cold Model", "Full Model"]
    colors = ["#f97306", "#107ab0", "#3f9b0b"]
    path = os.path.join(
        "figures", "site_experiment_scatters", f"cold_on_hot_scatter.svg"
    )
    with torch.no_grad():
        plot.fit_plot(
            models_list,
            cold_data,
            labels,
            path,
            unit=f"Q_LE [W/m²]",
            pointsize=6,
            alpha=0.2,
            colors=colors,
            show=True,
        )

    # Then, plot the hot, cold and balanced pipeline with hot data
    path = os.path.join(
        "figures", "site_experiment_scatters", f"hot_on_cold_scatter.svg"
    )
    with torch.no_grad():
        plot.fit_plot(
            models_list,
            hot_data,
            labels,
            path,
            unit=f"Q_LE [W/m²]",
            pointsize=6,
            alpha=0.2,
            colors=colors,
            show=True,
        )

    # Now: plot the dry, wet and balanced pipeline with wet data
    models_list = [dry_pipeline, wet_pipeline, balanced_pipeline]
    labels = ["Arid Model", "Humid Model", "Full Model"]
    colors = ["#e2ca76", "#047495", "#3f9b0b"]
    path = os.path.join(
        "figures", "site_experiment_scatters", f"wet_on_dry_scatter.svg"
    )
    with torch.no_grad():
        plot.fit_plot(
            models_list,
            wet_data,
            labels,
            path,
            unit=f"Q_LE [W/m²]",
            pointsize=6,
            alpha=0.2,
            colors=colors,
            show=True,
        )

    # Then, plot the dry, wet and balanced pipeline with dry data
    path = os.path.join(
        "figures", "site_experiment_scatters", f"dry_on_wet_scatter.svg"
    )
    with torch.no_grad():
        plot.fit_plot(
            models_list,
            dry_data,
            labels,
            path,
            unit=f"Q_LE [W/m²]",
            pointsize=6,
            alpha=0.2,
            colors=colors,
            show=True,
        )


def coeff_determination_ordered(test_datas, rs_model, seed=24):
    pipeline = pipelines.make_pipeline(rs_model, None, None)
    empirical_pipeline = pipelines.make_pipeline(None, None, None)
    gen = np.random.default_rng(seed)
    combined = []
    for site_data, site_name in test_datas:
        print(f"Evaluating site {site_name}:")
        # Take random sample
        test_indices = gen.choice(
            len(site_data),
            2000,
            replace=False,
        )
        site_data = np.array(site_data, dtype=object)[test_indices]

        our_coeff = eval.coefficient_of_determination(pipeline, site_data)
        empirical_coeff = eval.coefficient_of_determination(
            empirical_pipeline, site_data
        )

        print(f"Hybrid R² (ours) : {our_coeff}")
        print(f"Empirical R²: {empirical_coeff}")

        combined += site_data

    our_combined_coeff = eval.coefficient_of_determination(pipeline, combined)
    empirical_combined_coeff = eval.coefficient_of_determination(
        empirical_pipeline, combined
    )
    print(f"Evaluating Combined Data:")
    print(f"Hybrid R² (ours) : {our_combined_coeff}")
    print(f"Empirical R²: {empirical_combined_coeff}")


def coeff_determination(test_data, rs_model):
    pipeline = pipelines.make_pipeline(rs_model, None, None)
    empirical_pipeline = pipelines.make_pipeline(None, None, None)
    our_coeff = eval.coefficient_of_determination(pipeline, test_data)
    empirical_coeff = eval.coefficient_of_determination(empirical_pipeline, test_data)
    print(f"Hybrid R² (ours) : {our_coeff}")
    print(f"Empirical R²: {empirical_coeff}")


def evaluate_site_models(cold, hot, wet, dry, test_data, train_data):
    print(f"evaluating cold model on test data")
    coeff_determination_ordered(test_data, cold)
    print(f"evaluating hot model on test data")
    coeff_determination_ordered(test_data, hot)
    print(f"evaluating wet model on test data")
    coeff_determination_ordered(test_data, wet)
    print(f"evaluating dry model on test data")
    coeff_determination_ordered(test_data, dry)

    print(f"evaluating cold model on train data")
    coeff_determination_ordered(train_data, cold)
    print(f"evaluating hot model on train data")
    coeff_determination_ordered(train_data, hot)
    print(f"evaluating wet model on train data")
    coeff_determination_ordered(train_data, wet)
    print(f"evaluating dry model on train data")
    coeff_determination_ordered(train_data, dry)


# Objects and data needed for plotting
base_path = os.path.expanduser("~/epfl/semester_project/databases/T_C_PIPELINE_DATA/")
test_sites = ["AU-TTE", "CH-Cha", "CG-Tch", "FI-Sod"]
train_sites = ["CN-Du2", "US-Me2", "CH-Dav", "ES-Amo", "AU-Wac", "AU-Stp", "GF-Guy"]

best_rs_model = load_model("best_model")
hot_rs_model = load_model("hot_sites_only")
cold_rs_model = load_model("cold_sites_only")
dry_rs_model = load_model("dry_sites_only")
wet_rs_model = load_model("wet_sites_only")
# rs_model = load_model('default')

multiple_rs_models = [
    best_rs_model,
    load_model("best_model_2"),
    load_model("best_model_3"),
    load_model("best_model_4"),
]

# train_data, test_data = load_train_test_data('best_model')
# train_data, test_data = load_train_test_data('default')

# ordered_data = load_ordered_data('default')
# train_ordered_data = load_ordered_data_from_sites(train_sites, base_path)
# test_ordered_data = load_ordered_data_from_sites(test_sites, base_path)


# Call the desired functions generating plots
# generate_hp_tuning_loss_plot()
# generate_torch_viz_graph()
# generate_rs_slice_plots(rs_model)
# generate_Q_LE_slice_plots(rs_model)
# generate_scatter_plots(test_data, rs_model)

# interesting intervals by site:
# CN-Du2: 1500-1800
# US-Me2: 2400-2700?
# CH-Dav: very variable. maybe smtg like 6400-6700 to show different intensities
# ES-Amo: Anywhere is good, the empirical model so drastically overestimates
# Au-Wac: bad throughout, both the empirical and us underestimate, but we do a little more
# Au-Stp: Anywhere is good, same as ES-Amo but even more
# Gf-Guy: Most places, empirical underestimates, we do better.
# ranges = [(1500, 1800), (2400, 2700), (6400, 6700), (1500, 1800), (1500, 1800), (1500, 1800), (1500, 1800)]
# generate_Q_LE_timeseries(train_ordered_data, ranges, best_rs_model)
# test_ranges = [(1000, 1300), (1000, 1300), (1000, 1300), (1000, 1300)]
# generate_Q_LE_timeseries(test_ordered_data, test_ranges, best_rs_model)

# generate_site_experiment_plots(hot_rs_model, cold_rs_model, dry_rs_model, wet_rs_model, best_rs_model, test_ordered_data)

# generate_empirical_learning_scatter('best_model')
# generate_Q_LE_rs_sensitivity_plots()
# generate_multiple_model_plots(multiple_model_config_names)

# generate_scatter_plots(test_data, best_rs_model)

# generate_multiple_model_plots(multiple_rs_models)

# Call some evaluation functions
# coeff_determination_ordered(test_ordered_data, best_rs_model)
# evaluate_site_models(cold_rs_model, hot_rs_model, wet_rs_model, dry_rs_model, test_ordered_data, train_ordered_data)
# coeff_determination(test_data, best_rs_model)
