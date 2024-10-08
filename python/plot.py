import matplotlib.pyplot as plt
import numpy as np
import torch

def plot_losses(losses):
    plt.plot(losses, scaley="log")
    plt.show()

def plot_univariate_slice(model, f, x_0, idx, x_idx_range):
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
    plt.savefig(f"results/figures/univariate_plots/uv_idx_{idx}")
    plt.show()

    