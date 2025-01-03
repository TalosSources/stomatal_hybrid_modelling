import matplotlib.pyplot as plt
import numpy as np
import torch

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
            x_pos = int((0.7 + 0.3*i/len(coefficients))* len(smoothed_loss))  # Position slightly to the left of the rightmost point
            y_pos = smoothed_loss[x_pos] * (1.15 if i!=0 else 0.85)  # Position slightly above the curve
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
def fit_plot(model, f, points, model_name, show=False):
    # TODO: Make more beautiful, legible (add colors and shapes), 
    # to help distiguish several different scatters (add option for that),
    # and add many eloquent legend options 
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
    # TODO: like all others. more beautiful, more options, several timeseries,
    # legend options, other cool stuff
    # include: a way to do the 13. figure idea, namely a breakdown of the Q_LE time-serie among contexts.
    n = len(timeseries[0][0])
    r = range(n)
    for timeserie, subname in timeseries:
        plt.plot(r, timeserie, label=subname)

    plt.legend()

    plt.savefig(f"results/figures/time_serie_plots/{name}")
    if show:
        plt.show()

"""
Might be quite similar to a univariate_slice plot.
"""
def sensitivity_plot():
    # TODO
    ...
        
    