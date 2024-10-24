"""
In this file are differentiable functions corresponding to 
known relationships and formulas that can be used in a 
differentiable pipeline 
"""
import torch

def Vmax():
    # TODO: Should it be done?
    ...

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