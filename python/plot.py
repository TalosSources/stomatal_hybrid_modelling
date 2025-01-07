import matplotlib.pyplot as plt
import numpy as np
import torch

from scipy.stats import multivariate_normal

import utils

# TODO: make more beautiful, give more options (legend etc.), make shaded min-high values, smoothing?
# or perhaps just use wandb for this one?
def plot_losses(losses, labels, coefficients=None, minmax=True, smoothing=0.5):
    plt.close('all')

    """
    Plot one or several given series of losses (losses),
    each using a different color. 
    If minmax is set, show the min and max values of the loss as a shaded area 
    around the main plot lines, of the corresponding color.
    Label the plots using labels, and show an appropriate legend.

    Parameters:
    - losses: List of lists or arrays, where each entry is a loss series.
    - labels: List of labels corresponding to each loss series.
    - minmax: Boolean, whether to show min/max shaded areas. Default is True.
    - smoothing: Float between 0 and 1 for exponential smoothing. Default is 0.5.
    """
    def smooth_data(data, alpha):
        smoothed = []
        last = data[0]  # Initialize with the first value
        for point in data:
            smoothed_val = alpha * point + (1 - alpha) * last
            smoothed.append(smoothed_val)
            last = smoothed_val
        return smoothed

    plt.figure(figsize=(10, 6))

    for i, loss in enumerate(losses):
        label = labels[i] if i < len(labels) else f"Series {i + 1}"

        if smoothing:
            smoothed_loss = smooth_data(loss, smoothing)
        else:
            smoothed_loss = loss

        line, = plt.plot(smoothed_loss, label=label, linewidth=2)
        color = line.get_color()

        if coefficients is not None and i < len(coefficients):
            coefficient = coefficients[i]
            x_pos = int((0.69 + (0. if i==0 else 0.5)/len(coefficients))* len(smoothed_loss))  # Position slightly to the left of the rightmost point
            y_pos = smoothed_loss[x_pos] * (1.2 if i==2 else (0.81 if i==1 else 0.80))  # Position slightly above the curve
            plt.text(x_pos, y_pos, f"R² = {coefficient:.2f}",
                     color=color, fontsize=12, fontweight='bold', verticalalignment='center')

        if minmax:
            min_vals = np.minimum.accumulate(loss)
            max_vals = np.maximum.accumulate(loss)
            plt.fill_between(range(len(loss)), min_vals, max_vals, alpha=0.2, color=color)

    plt.yscale("log")
    plt.xlabel("Epochs")
    plt.ylabel("Loss (Log Scale)")
    plt.title("Loss Over Time")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    #plt.tight_layout()
    plt.show()

def plot_univariate_slices_subplots(models, x_0s, x_0_labels, idxes, steps, model_labels, res_label, path, x_0_suppl=False, show=False):
    plt.close('all')
    
    n_idxes = len(idxes)
    n_x_0s = len(x_0s)
    fig, axs = plt.subplots(n_idxes, n_x_0s, figsize=(4 * n_x_0s, 3 * n_idxes), squeeze=False)

    for i, (idx, idx_range, idx_label) in enumerate(idxes):
        # create the values we explore for dimension idx
        range_start, range_stop = idx_range
        step = (range_stop - range_start) / steps
        x_range = torch.arange(range_start, range_stop, step)

        for j, (x_0, x_0_label) in enumerate(zip(x_0s, x_0_labels)):
            ax = axs[i, j]
            # create an array of copies of each x_0, where we replace idx's value by one of the explored values
            if x_0_suppl:
                x_0, datapoint = x_0
            x_inputs = torch.tile(x_0, (steps, 1))
            x_inputs[:, idx] = x_range[:]

            for model, label in zip(models, model_labels):
                y_range = model(x_inputs, datapoint) if x_0_suppl else model(x_inputs)
                label = label if (x_0_label == "") else f"{label}"
                ax.plot(x_range, y_range, label=label, linewidth=2)

            #ax.set_title(f"{idx_label}", fontsize=10)  # Dimension label for each subplot
            ax.set_xlabel(idx_label, fontsize=8)
            ax.set_ylabel(res_label, fontsize=8)
            ax.grid(True, linestyle='--', alpha=0.6)
            ax.set_box_aspect(1)

    for j, x_0_label in enumerate(x_0_labels):
        fig.text(
            0.5 / n_x_0s + j / n_x_0s,  # x-coordinate (centered above each column)
            0.98,  # y-coordinate (just above the subplots)
            x_0_label,
            ha='center',
            fontsize=10,
        )

    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', fontsize=10)
    plt.tight_layout(rect=[0, 0, 1, 1])  # Leave space for legend
    
    #plt.ylabel(res_label)

    plt.savefig(path)

    if show:
        plt.show()

def plot_univariate_slices(models, x_0s, x_0_labels, idx, idx_label, range, steps, model_labels, res_label, path, show=False):
    plt.close('all')
    plt.figure(figsize=(10, 6))

    # create the values we explore for dimension idx
    range_start, range_stop = range
    step = (range_stop - range_start) / steps
    x_range = torch.arange(range_start, range_stop, step)

    for x_0, x_0_label in zip(x_0s, x_0_labels):
        # create an array of copies of each x_0, where we replace idx's value by one of the explored values
        x_inputs = torch.tile(x_0, (steps, 1))
        x_inputs[:, idx] = x_range[:]

        for model, label in zip(models, model_labels):
            y_range = model(x_inputs)
            label = label if (x_0_label == "") else f"{label} ({x_0_label})"
            plt.plot(x_range, y_range, label=label, linewidth=2)

    plt.xlabel(idx_label)
    plt.ylabel(res_label)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    #plt.title("Univariate slice of rs models")

    plt.savefig(path)

    if show:
        plt.show()

def plot_univariate_slice(model, f, x_0, idx, x_idx_range, model_name, show=False):
    # TODO: Make a much more complete and general version. 
    # Make it beautiful, and potentially containing more information, several slices or several sites/timesteps? (=x0) 
    # (take inspiration from paper)
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
def fit_plot(models, points, labels, path, unit, obs_label='Observed', pred_label='Predicted', plot_range=None, pointsize=10, alpha=0.3, colors=None, show=False):
    # TODO: Make more beautiful, legible (add colors and shapes), 
    # to help distiguish several different scatters (add option for that),
    # and add many eloquent legend options 
    plt.close('all')
    plt.figure(figsize=[7, 5])

    if colors is None:
        colors = ['#1f77b4']*len(models)

    # TODO: crop to valid region, > 500 is broken flux anyway?

    use_grouped = False

    if use_grouped:
        grouped, _ =  utils.group_by_ctx_id(points)
        x_points = []
        for group in grouped:
            x_points = x_points + np.array([y.numpy() for (_, y) in group])
    else:
        x_points = np.array([y.numpy() for (_,y) in points])
    # Compute the points for the plot
    for model, label, color in zip(models, labels, colors):
        if use_grouped:
            y_points = []
            for group in grouped:
                data = utils.make_batch(group)
                outs = model(data).numpy()
                y_points = y_points + list(outs)
            y_points = np.array(y_points).flatten()
        else:
            y_points = np.array([model(x).numpy() for (x,_) in points]).flatten()

        scatter = plt.scatter(x_points, y_points, label=label, s=pointsize, alpha=0.6, color=color)

        # Compute the 2D Gaussian fit
        data = np.vstack((x_points, y_points)).T
        mean = np.mean(data, axis=0)
        cov = np.cov(data, rowvar=False)
        gaussian = multivariate_normal(mean=mean, cov=cov, allow_singular=True)

        # Create grid for contour plot
        min_value = min(x_points.min(), y_points.min())
        max_value = max(x_points.max(), y_points.max())

        x, y = np.meshgrid(
            np.linspace(min_value, max_value, 100),
            np.linspace(min_value, max_value, 100)
        )
        pos = np.dstack((x, y))
        z = gaussian.pdf(pos)

        # Plot filled contour for a specific density threshold
        plt.contourf(x, y, z, levels=[z.max() * 0.1, z.max()], colors=[color], alpha=0.3)

        # Plot contour for a specific density threshold
        plt.contour(x, y, z, levels=[z.max() * 0.1], colors=[color], linestyles='solid', linewidth=1)

    # Plot the identity line x = y
    max_value = max(max(x_points), max(y_points))
    min_value = min(min(x_points), min(y_points))
    plt.plot([min_value, max_value], [min_value, max_value], 'r--', label="x = y (Ideal fit)")

    # Set axis limits and make sure both axes start from 0 with the same scale
    min_plot = 0
    max_plot = max_value
    if plot_range is not None:
        min_plot, max_plot = plot_range
    plt.xlim(min_plot, max_plot)
    plt.ylim(min_plot, max_plot)
    plt.gca().set_aspect('equal', adjustable='box')

    # Add labels and title
    plt.xlabel(f"{obs_label} {unit}")
    plt.ylabel(f"{pred_label} {unit}")
    plt.legend()

    plt.tight_layout()

    plt.savefig(path, dpi=300) # TODO: dpi=300 and also on all other figures
    if show:
        plt.show()

def time_series_plot(timeseries, labels, path, linestyles, widths=None, alphas=None, x_label=None, y_label=None, ranges=None, show=False):
    # TODO: like all others. more beautiful, more options, several timeseries,
    # legend options, other cool stuff
    # include: a way to do the 13. figure idea, namely a breakdown of the Q_LE time-serie among contexts.
    plt.close('all')
    plt.figure(figsize=(10, 4))
    n = len(timeseries[0])
    if ranges is None:
        ranges = [range(n)] * len(timeseries)
    if widths is None:
        widths = [1] * len(timeseries)
    if alphas is None:
        alphas = [0.8] * len(timeseries)
    for timeserie, r, label, linestyle, width, alpha in zip(timeseries, ranges, labels, linestyles, widths, alphas):
        plt.plot(r, timeserie, label=label, linestyle=linestyle, linewidth=width, alpha=alpha)

    if x_label is not None:   
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)
    plt.legend()
    plt.tight_layout

    plt.savefig(path)
    if show:
        plt.show()

"""
Might be quite similar to a univariate_slice plot.
"""
def sensitivity_plot():
    # TODO
    ...
        
    