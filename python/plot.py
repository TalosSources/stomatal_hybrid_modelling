import matplotlib.pyplot as plt
import numpy as np
import torch

def plot_losses(losses):
    plt.close('all')

    plt.plot(losses, scaley="log")
    plt.show()

def plot_univariate_slice(model, f, x_0, idx, x_idx_range, model_name, show=False):
    plt.close('all')

    y_range_f = np.zeros_like(x_idx_range)
    y_range_model = np.zeros_like(x_idx_range)
    x = torch.clone(x_0)
    for i in range(y_range_f.shape[0]):
        x[idx] = x_idx_range[i]
        y_range_f[i] = f(x)
        y_range_model[i] = model(x)

    print(f"x_idx_range={x_idx_range}")
    print(f"y_f={y_range_f}")
    print(f"y_model={y_range_model}")

    plt.plot(x_idx_range, y_range_f, label="f plot")
    plt.plot(x_idx_range, y_range_model, label="model plot")
    plt.savefig(f"results/figures/univariate_plots/{model_name}_idx_{idx}")
    
    if show:
        plt.show()

"""
model and f are real functions, and we plot the fit of the model to f.
For each point in points, we plot the 2D point (f(point), model(point)).
If the functions evaluate to the same values for a point, it should lie
on the x = y line. 
"""
def fit_plot(model, f, points, model_name, show=False):
    plt.close('all')

    # Compute the points for the plot
    x_points = [f(x).detach().numpy() for x in points]
    y_points = [model(x).detach().numpy() for x in points]
    
    plt.scatter(x_points, y_points, label="Fit plot")

    # Plot the identity line x = y
    max_value = max(max(x_points), max(y_points))
    min_value = min(min(x_points), min(y_points))
    plt.plot([min_value, max_value], [min_value, max_value], 'r--', label="x = y (Ideal fit)")

    # Set axis limits and make sure both axes start from 0 with the same scale
    plt.xlim(0, max_value)
    plt.ylim(0, max_value)
    plt.gca().set_aspect('equal', adjustable='box')

    # Add labels and title
    plt.xlabel("f(x)")
    plt.ylabel("model(x)")
    plt.title(f"Fit plot for {model_name}")
    plt.legend()

    plt.savefig(f"results/figures/fit_plots/{model_name}_{len(points)}_points")
    if show:
        plt.show()

def time_series_plot(timeseries, name="timeseries", show=False):
    n = len(timeseries[0][0])
    r = range(n)
    for timeserie, subname in timeseries:
        #print(f"plotting {subname}... (shape={len(timeserie)})")
        #print(f"timeserie: {timeserie}")
        plt.plot(r, timeserie, label=subname)

    plt.legend()

    plt.savefig(f"results/figures/time_serie_plots/{name}")
    if show:
        plt.show()
        
    