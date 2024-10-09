* what is written in comments around gsco2
* find where photosynthesis_biochemical is being called, the module just before and the module just after, and translate both of them as pytorch modules.global datasets of hourly 
[
    In CO2_Concentration.m: actually it's just a wrapper of pb.
    In Canopy_Resistence_An_Evolution.m: fzero is called just before pb both times.
    it's not really calling any function before or after actually, it's just calling a wrapper of itself.
    let's translate CO2_Concentration.m then, but I don't know an equivalent to fzero right now.
]
* look for global datasets of stomatal conductance. possible starting point: the saved paper about optimal stomatal behaviour
* go through the "global datasets of hourly... " paper and identify what they do with VCmax. find other parameters to tune in the same way they tune VCMax (using neural nets, random forests, ...). in the paper: they construct a model for the slope a1 of gsco2, which we could do as well. we may tune Vm or Vmax. actually Vmax and Vcmax may be the same. we might remove the line 
* -> create a NN for VMax, train both at the same time

* if time: take the ... paper's model for evapotranspiration, but replace gsco2 by our gsco2, use their equation to obtain evapotranspiration. rs = stomatal resistance (or conductance). other ones I can assume arbitrary values? yes, and also for Q_LE.
* try to train a stochastic model? bayesian learning?

can contact Akash until friday morning, and then not until 08.10


# TODO Starting 08.10
* make a fit plot function, smtg smtg coefficient of determination, smtg smtg R²
* Setup MATLAB
* check that I can learn the VCMax empirical function
* find where ra is coming from in pb.m
* write the logic in torch for going from gsco2 to Q_LE and Q_H
* obtain Q_LE (and later Q_H) dataset [PENDING ON DATA]
* train a hybrid model of GSCO2, and VCMax (large pipeline) using a loss function including Q_LE like in the evapotranspiration dataset  [PENDING ON DATA]
* understand how they constrain the hybrid models, see if I can do it as well
* write report (?)
* Explore Bayesian Deep Learning (links saved in semester_project)
