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