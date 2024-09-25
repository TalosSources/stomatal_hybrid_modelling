* what is written in comments around gsco2
* find where photosynthesis_biochemical is being called, the module just before and the module just after, and translate both of them as pytorch modules.global datasets of hourly
* look for global datasets of stomatal conductance. possible starting point: the saved paper about optimal stomatal behaviour
* go through the "global datasets of hourly... " paper and identify what they do with VCmax. find other parameters to tune in the same way they tune VCMax (using neural nets, random forests, ...). in the paper: they construct a model for the slope a1 of gsco2, which we could do as well. we may tune Vm or Vmax. actually Vmax and Vcmax may be the same. we might remove the line 
* -> create a NN for VMax, train both at the same time

* if time: take the ... paper's model for evapotranspiration, but replace gsco2 by our gsco2, use their equation to obtain evapotranspiration. rs = stomatal resistance (or conductance). other ones I can assume arbitrary values? yes, and also for Q_LE.
* try to train a stochastic model? bayesian learning?

can contact Akash until friday morning, and then not until 08.10