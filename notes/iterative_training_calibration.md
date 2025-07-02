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
~* add a better site support. Need to answer some questions:
    1. do we train on all sites? a selection?
    2. do we predict on all sites? a selection? makes sense to predict only on sites we use to train (otherwise we see no progress?)
    3. do we evaluate on all sites?
    In fact, we need, when we train in python, to use sites that were predicted using the most recent version of rs_model. Otherwise, we're not making the most progress we can. Then, to evaluate the progress, we need to evaluate on the same as well.
    In fact, to hope to guarantee any progression, we need to fix the site list over the whole iterative training procedure.~
* perhaps make the script more parametrizable-ed. For example, the rs model path in matlab shouldn't ideally be hard-coded

## Then
* hyperparams again?

## Goal
* do the above (technical), and check everything works
* ensure MAE is going down (otherwise a bit strange...)
* select sites
* In meeting, present the running pipeline. Find a way to run it (cluster? on my own computer? but takes a while. GPU?)
* ask what results we want (what plots, what metrics, etc.). What do we have more using T&C? new variables can be plotted? or just more accurate results?
* run the iterative training calibration
* make some technical report with the outputs
* clean the code and repo for Son
    * refactor -> change/fix the t_and_c path, make sure everything is a parameter
    * write README, check requirements.txt is complete

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