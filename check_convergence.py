import os
import sys
import scipy
import numpy as np
import csv

# NOTE: Should we check the predicted Q_LE values by T&C?
# or maybe load the previous and current rs_model, and simply test for some slice-like values
# -> test for An convergence. consider MAE, to avoid peeks. relative error
prev_tandc_output_path = sys.argv[2]
current_tandc_output_path = sys.argv[1]
sites_csv_path = sys.argv[3]

threshold = 1e-5  # NOTE: How to choose it?

losses = []

print(f"received args={sys.argv}", file=sys.stderr)

# for test_site in test_sites:
# result_file = f"Results_{test_site}.mat"
# prev_path = os.path.join(prev_tandc_output_path, result_file)
# current_path = os.path.join(current_tandc_output_path, result_file)


sites = []
with open(os.path.expanduser(sites_csv_path), newline="") as sites_csv_file:
    reader = csv.reader(sites_csv_file)
    first = True
    for line in reader:
        if first:
            first = False
            continue
        else:
            sites.append(line[0])

for site in sites:
    prev_path = os.path.join(prev_tandc_output_path, f"Results_{site}.mat")
    current_path = os.path.join(current_tandc_output_path, f"Results_{site}.mat")

    prev_data = scipy.io.loadmat(prev_path)
    current_data = scipy.io.loadmat(current_path)

    prev_An_L = prev_data["An_L"]
    prev_An_H = prev_data["An_H"]
    current_An_L = current_data["An_L"]
    current_An_H = current_data["An_H"]

    loss_H = np.abs(prev_An_H - current_An_H).mean()
    loss_L = np.abs(prev_An_L - current_An_L).mean()
    losses.append(loss_H)
    losses.append(loss_L)

mean_loss = np.array(losses).mean()
print(f"CHECKING CONVERGENCE: Found 'An' loss = {mean_loss}", file=sys.stderr)
if mean_loss > threshold:
    print("false", sep="", end="", file=sys.stdout)
else:
    print("true", sep="", end="", file=sys.stdout)
