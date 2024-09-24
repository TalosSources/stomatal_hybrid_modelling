
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
### Measurements
### Physical constants
### Empirical parameters (to tune)



# gsCO2 dependancy
(in the end, clamping with go)
* gsCO2(go, a1, An, Pre, Cc, GAM, Ds, Do)