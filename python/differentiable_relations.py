"""
In this file are differentiable functions corresponding to 
known relationships and formulas that can be used in a 
differentiable pipeline 
"""
import numpy as np

def Vmax():
    # NOTE: Should it be done?
    ...


# Constants
cp = 1004.834
specific_heat_air = 1.013e3  # J/kg/°C
air_pressure_at_sea_level = 1.013e5  # Pa
specific_gravity_water_vapor = 0.622  # unitless
b0 = 1.91846e6
b1 = -33.91
zero_celsius_in_kelvin = 273.15   

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

Notes:
* For now, we assume Tmean is 15 but it is to be replaced by a predictor
"""
def Q_LE(rs, ra, Rn, QG, Ds, Ts, Tmean=15):
    lambda_ = b0 * np.power(Ts / (Ts + b1), 2)
    gamma = specific_heat_air * air_pressure_at_sea_level / (specific_gravity_water_vapor * lambda_)
    roa = air_pressure_at_sea_level / (287.05 * (Ts + zero_celsius_in_kelvin))
    sc = compute_sc(Tmean)
    return (sc * (Rn - QG) + roa*cp*Ds / ra ) / (sc + gamma*(1 + rs/ra))

def compute_rs(Q_LE, ra, sc, Rn, QG, roa, cp, es, ea, gamma):
    return (ra*sc*(Rn - QG) + roa*cp*(es -ea) - ra*Q_LE*(sc + gamma) ) / ( gamma * Q_LE )

def compute_rs_2(gsCO2, Tf, Ts, Pre, Pre0):
    rsCO2=1/gsCO2
    rsH20 = (rsCO2/1.64)*(1e6)
    return rsH20*(Tf*Pre)/(0.0224*(Ts+273.15)*Pre0)

# Taken from https://edis.ifas.ufl.edu/publication/AE459
def compute_sc(Tmean):
    return 4098 * (0.6108 * np.exp(17.27*Tmean / (Tmean + 237.3))) / np.power(Tmean + 273.3, 2) # NOTE: remplace np by torch?


#def Q_LE(gsCO2, Tf, Ts, Pre, Pre0, ra, sc, Rn, QG, roa, cp, es, ea, gamma):
#    ... # NOTE: Other version with other dependancies, probably not useful
