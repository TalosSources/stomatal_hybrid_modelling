"""
In this file are differentiable functions corresponding to 
known relationships and formulas that can be used in a 
differentiable pipeline 
"""
import torch

def Vmax():
    # TODO: Should it be done?
    ...

"""
Parameters:
* ra is given (but should it? in the paper it's predicted) [Mostly OK]
* rs is a pb output (but also a predictor?) [OK]
* sc ??? slope of the saturation vapor pressure-temperature relationship (kPaC^-1) [ASK]
* Rn appears in Observations, and PredictorKeys??? called Net Radiation [Probably OK, Ask confirmation]
* QG ??? called Soil Heat Flux (same unit as QLE), may map to G (predictors and observations)? [Ask confirmation]
* roa: mean air density at constant pressure (kg*m^-3), perhaps ros in Predictors? [Ask]
* cp: specific heat of dry air at constant pressure = 1004.834 (Jkg^(-1) C^(-1)) [OK]
* es "esat" appears in Observations (and PredictorKeys???) [Ask clarification]
* ea appears in Observations (and PredictorKeys???)
-> es-ea: vapor pressure deficit (VPD) of air (kPa)
-> actually, probably Ds = (es - ea), and Ds is a predictor
* gamma: psychrometric constant (kPaC^-1) [ASK]
"""
def Q_LE(rs, ra, sc, Rn, QG, roa, cp, es, ea, gamma):
    return ( sc * (Rn - QG) + roa*cp*(es - ea) / ra ) / (sc + gamma*(1 + rs/ra))

def rs(Q_LE, ra, sc, Rn, QG, roa, cp, es, ea, gamma):
    return (ra*sc*(Rn - QG) + roa*cp*(es -ea) - ra*Q_LE*(sc + gamma) ) / ( gamma * Q_LE )

def rs(gsCO2, Tf, Ts, Pre, Pre0):
    rsCO2=1/gsCO2
    rsH20 = (rsCO2/1.64)*(1e6)
    return rsH20*(Tf*Pre)/(0.0224*(Ts+273.15)*Pre0)

def Q_LE(gsCO2, Tf, Ts, Pre, Pre0, ra, sc, Rn, QG, roa, cp, es, ea, gamma):
    ... # TODO

def Q_LE(gsCO2):
    return gsCO2 # TODO: Remove. simply for error-free run