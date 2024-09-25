
flag CT selects a technique

# Flow
* define Pre0, Pre, Tf, IPAR, Tref, R (constants)
* prepare Ta, IPAR, Cc, Oa, Csl, go, Ts_k
* compute ra(Ta, Pre0, Tf, Pre), rb(Ta, Pre0, Tf, Pre), rmes(gmes) using some midly complex formula
* select between 2 techniques to compute:
    constants Hd, Ha, DS (they have a very strange control flow)
    complex kT(Ha, Ts_k, Tref, R, DS, Hd)
    simple Vm(Vmax, kT)
    simple Jmax(Vmax, rjv)
    simple Jm(Jmax, kT)
and so on and so on

## Parameters
We must first decide what parameters can be modified (it makes no sense to modify the water fusion temperature of water in Kelvins, or a measurement of the carbon density in the atmosphere), and then it would be good to have validity range for them. Perhaps they should be positive, larger or smaller than 1, or something akin to that.
### Measurements
### Physical constants
### Empirical parameters (to tune)
candidates: activation energy (Ha), entropy factors (DS)
magic numbers: (4.57, 0.0224)
misc: (FI, Oa)
also: some formulas seem arbitrary / semi-empirical. I guess they may be replaced by simple parametrized models (small neural networks, adding some tunable parameters). for example: all the mess starting at ANS_TEMP



# gsCO2 dependancy
(in the end, clamping with go)
* gsCO2(go, a1, An, Pre, Cc, GAM, Ds, Do)

# Plan
We assume we get ground truth for gsCO2 for given arguments/measurements of the function. We need to need to create a list of tunable parameters in the pipeline, feed it to an optimizer, then repeatedly call the function with the data, use a simple loss such as MSE, and 