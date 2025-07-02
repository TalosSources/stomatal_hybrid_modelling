To correct for a bias in the training

# Plan
* train the model as usual using data from vanilla T&C
* repeat until convergence:
    * plug the trained model in T&C
    * generate new data by running T&C on all training sites
    * train the model using this new data

Can this be automated?
With smart handling of model and data paths, a bash script could easily iteratively start training and prediction rounds.
The only tricky part would be checking convergence. This would require saving the new data/model elsewhere than the old, so that both can co-exist and a comparison can be made. Or the comparison can be made before saving in memory. In any case, it would make sense to check only a small subset of the entire data for convergence, which would be equivalent to checking everything if the subset is representative and would save a lot of time/space. 
And to perform the actual check, it would be more convenient to do it in python, and output a flag that's readable by bash.

## Convergence
Can be defined in terms of model weights convergence, or predictions convergence. It seems the former always has significant stochasticity even when training goes well, so latter might make more sense.

# Roadmap
* remove prints, make the training understandable
* perhaps show some "convergence error" values, even a plot?

## Then
* hyperparams again?

## Goal
* do the above (technical), and check everything works
* ensure MAE is going down (otherwise a bit strange...)
* select sites
* clean the code and repo for Son
    * an inference script? for T&C, that's quite automatic
    * refactor -> change/fix the t_and_c path, make sure everything is a parameter
    * write README, check requirements.txt is complete
    * readme should be usable by Sara and Akash as well, so that they can make modifications, change sites etc., make plots and analysis.
* check the installation and script works for a fresh system 
* run the iterative training calibration with the chosen sites

# Perf
## Knobs
* python: data size, epochs, batch size, network size
* matlab: iters, perhaps also data size
* sites count
## modified
* removed the python part
* changed iters and data size in python
* changed iters in matlab

# Problems
* no convergence -> consider using r_s predictions instead of An, or something else

# General notes
the script assumes that T&C was run before, with the modified final parameters. Perhaps, that should be mentioned somewhere, or the script should (offer the option to) run it before running python