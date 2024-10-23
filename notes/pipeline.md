# Description of the differentiable pipeline
* We have access to ground truth for Q_LE and another one
* We have access to some physical measurements that we treat as our predictors
* We know gsCO2 and ra, which are hard to measure directly, are somehow a function of the predictors 
* We know of simple and differentiable functions from gsCO2 and ra to Q_LE and the other

With all that:
* we create a regressor model from the predictors to gsCO2 and ra
* For each predictor point set, we predict gsCO2 and ra, then differentiably obtain Q_LE, and use some error loss to train the regressor. All the operations downstream of the regressor must be differentiable.
 