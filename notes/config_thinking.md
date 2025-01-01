# Ablation study
## Methodology
some parameters are interesting, others less.
Notably, the interesting parameters relate 

the plan can be:
For each interesting set of parameters, perform a hp-tuning run.
Then, use the best hps to perform a large run.

## Study
### Hyperparameters
* lr
* weight decay
* model params:
    * nHidden [2,3,4,5]?
    * hidden_size [32, 64]?
    * batch normalization / not
* ev (but probably not): batch size

### Interesting Parameters
* what variable we model: [rs, gsCO2]
* potentially: what sites to use, so that we can say: 
    using dry sites, performs poorly on wet, vice-versa, but using variety of sites -> performs well on unknown sites
* whether to exclude datapoints with Q_LE < 0
* eventually: whether to use softplus (to show it's bad without it)

It will be complicated to test all possible combinations, both for hps, and for interesting params. I have to make a distinction between params for which I want to find the optimal values, and params for which I just want to make a point for the report. For the latter, I simply test possible values while the other params are fixed, and I report the results.
Actually, let's assume that for each of these choices, one is better than the other. Then, I need to run this: one run with all the best params, which may actually be the same as the best_model run (or not, for reasons such as 'number of epochs' / 'sample points' / 'I found better tweaks in the meantime and want the best result possible in the end'), and then one run for each of the sub-optimal interesting param choice (with all other params being their optimal version). With the current version, that's 1+sum[n_options - 1] = 5 I think? And for each run, we do hp tuning, which requires running about 10-20 cross-validation, which each require k training runs. If we use a lot of data-points, a training run could require 1000 epochs... putting everything together, with k=5, that'd require about 500 times this 1000 epoch run... which seems a bit too much, honestly. I think, I need to lower this number. 
Measures to lower this number:
1. combine nHidden and hidden_size: test specific tuples of those, not all possible tuples.
2. Don't perform the hp-tuning for each experiment. HP-tuning is for me, not the report. I'm using it to get good hps and obtain the best performance I can, I will briefly mention it in the report, but probably won't provide plots/details on this, or minimally so, as this is not the topic of my project. The most sensible thing might be to do HP-tuning on a "representative" run, using what I suspect might be the optimal param set, then run the experiments with this HP-set. Then, I might perform another, potentially expensive hp-tuning for my final all-out experiment. -> this means I don't need to annoyingly combine hp-tuning and training.

Note: If I'm not sure for all interesting params, which is the best, I should run in priority params for which I'm least confident, so that I can observer the results, and adjust its value as the optimal one when I run the other params.