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
* make a fit plot function, smtg smtg coefficient of determination, smtg smtg R² [DONE]
* Setup MATLAB [DONE]
* check that I can learn the VCMax empirical function [Sort of DONE]
* find where ra is coming from in pb.m [Probably DONE]
* write the logic in torch for going from gsco2 to Q_LE and Q_H
* obtain Q_LE (and later Q_H) dataset [PENDING ON DATA]
* train a hybrid model of GSCO2, and VCMax (large pipeline) using a loss function including Q_LE like in the evapotranspiration dataset  [PENDING ON DATA]
* understand how they constrain the hybrid models, see if I can do it as well
* write report
* Explore Bayesian Deep Learning (links saved in semester_project)

# TODO Starting 16.10
* start the gsCO2 -> Q_LE pipeline, use the Sites in CH data to train an actual gsCO2 model, evaluate it.
* find predictors for VCMax, check whether we have them already or we need external data, in which case ask Akash. Check if the paper (hourly...) has reference data. Maybe check other papers related to VCMax, come up with a set of predictors to include. [LOOK AT A FEW PAPERS; MAKE A LIST OF PREDICTORS WITH THEIR FULL NAMES; ASK AKASH]
* Model VMax as a network directly dependant on the above predictors, replace it in pb.py module.
* train the gsCO2 model for each site sequentially. Train all the data for the same site sequentially, then move to another site, instead of mixing data. [Compare methods? sequentially vs not-sequentially]. [important: some state variable depend on previous model steps, for example soil moisture today depends on soil moisture yesterday, so in practice the model is always run sequentially. Think about that, does it make sense to shuffle the data in this context? research markovian networks? maybe recurrent models [READ ABOUT THEM]? so also test including previous state data in the input and shuffle, compare results]
* [If Time] train the VMax and gsCO2 at the same time, evaluate the whole pipeline

# TODO Starting 06.11
* some stuff in the data file
* use gsco2 in the model files to check the pipeline has correct magnitude

# Later
* learn for different sites (other biomes), different sequence in which sites are trained

# 13.11 Notes
* idea: average over one day (48 time-steps in the data array) to reduce the noise
* give 2 time series plots: predicted vs observed LE
* ask akash to inspect my code (give him indications)
* see if LE_CORR improves things
* ask akash for dry sites (US, ...)


# 22.11 Notes
* test the inverted PM equation rs(Q_LE, ra, ...), to see if the obtained values match the given rs.
* come up with smart ways to pinpoint where the mistake in the pipeline is
* test that the PM function and the inverted PM are reciprocal (f^-1 (f(x)))

# 06.12 Notes
* For Matching Q_LE: it's a sum of sunH, sunL, shdH, shdL values, and also other things we don't compute, so we need to take all that into account
* We are calculating T_sunH and we're comparing with the sum of above values.
* We need to run our pipeline for shdH, we get T_shdH, T_H = T_sunH + T_shdH. (would be the Q_LE but converted from W/m² to mm/H (mm of water, look online?))
* Next, we compute ET =  sum(T_H+EIn_H,2) + sum(T_L+EIn_L,2) +  EG +  ELitter + ESN + ESN_In + EWAT +  EICE+ EIn_urb + EIn_rock ;  %% [mm/h]
* all the ones we need (all execpt T_H) should be in the predictors
* I could also use the pipeline to obtain T_L (that would amount to running it 4 times)

## Divergences
IPAR: "IPAR":"PAR_sun_H_final" is wrong
Psi_L: suspect -> corresponds to the previous one
Vmax is slightly 


# 11.12 Notes
* Aggregate values of one day to compute the loss function (to smooth out the noise)
* Use the correct Vmax (not a constant) ASK AKASH
* then fix psi_l etc.

* make sure there's some sensitivity in training: changing the predicted gsco2 does change the final pipeline output?
* perhaps plot it (the above)
* think about what to show in the report: what figures. (scatter plot, violin plot, compare pure T&C with hybrid, loss plot (loss goes down), hyperparameter tuning (different learning rates, different models), ablation study (test to remove some of the choices we made, e.g. to train with daily averages))

# 22.12 Notes
## Questions: 
* Ask about the time-step precise info: when does it start in the day/year?
## Notes
* Think about invalid data: If some given rs are inf, or other values are inf or nan, consider removing them
* Dates appear in results['Date'], a big array containing dates for each timestep. I need to understand and parse the encoding, can be useful to provide analysis about variations over the course of a day/year.
* need to ping Akash by email for notification
* consider accelerating with CUDA
* IDEA TO TRY: first, train the model with the empirical function, as a baseline. Forces it to be somewhat coherent. Then, train from it with the pipeline.
* weight decay
* layer normalization?
* write evaluation functions. It seem if I can learn anything, it will be hard to spot just by staring at outputs/lossPlots. I need to make evaluations before/after training, and against empirical models.
* Do is a constant. Does it need to be in the inputs of pb? doesn't it just makes computations more expensive?

## Matlab Code for Date:
The invert of one of those procedures should give me the human-readable dates from the data numbers.
### First Option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  JULIAN_DAY              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[jDay]= JULIAN_DAY(D)
%%%% INPUT
%%% Datam %% [Yr, MO, DA, HR]
% Determine the julian day of the current time
days = [31 28 31 30 31 30 31 31 30 31 30 31];
nowYR=D(1); nowMO=D(2); nowDA=D(3); %nowHR =D(4);
if(nowMO==1)
    jDay = nowDA;
elseif(nowMO==2)
    jDay = days(1) + nowDA;
else
    jDay = sum(days(1:(nowMO-1))) + nowDA;
    if(mod(nowYR,4)==0)
        if(mod(nowYR,400)==0)
            jDay = jDay + 1;
        elseif(mod(nowYR,100)~=0)
            jDay = jDay + 1;
        end
    end
end
return

### Second option
%% Extract and process the dates
date_start=datenum(num2str(data(1,2)),'yyyymmddHHMM');
date_end=datenum(num2str(data(end,2)),'yyyymmddHHMM');
% Hourly Date
Date = date_start:1/24:date_end;
Date=Date';

