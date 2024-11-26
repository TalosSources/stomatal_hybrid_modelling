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
air_pressure_at_sea_level_kPa = air_pressure_at_sea_level / 1000 # kPa
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
* I think the given Ds is in Pa, while it should be in kPa here
* this formula computes LE in W/m², which seems to be the ground truth values from T&C
* ra and rs are probably of the right unit (except if stuff is passed by reference and modified)
--> actually tensors seem to be partly passed by reference. check ra, rs and stuff keep their values
* Rn and QG are probably okay
* Ts is the important part:
    * sc should be kPa/C (is C for celsius??)
    * roa should be kg/m³
"""
def Q_LE(rs, ra, Rn, QG, Ds, Ts):

    sc_kPa, roa, Ds_kPa, gamma = compute_missing_variables(Ts, Ds)
    
    q_le = (sc_kPa * (Rn - QG) + roa*cp*Ds_kPa / ra ) / (sc_kPa + gamma*(1 + rs/ra)) # flux, W/m²
    #print(f"obtaining Q_LE={q_le}")
    return q_le

def compute_rs(Q_LE, ra, Ts, Ds, Rn, QG):
    sc_kPa, roa, Ds_kPa, gamma = compute_missing_variables(Ts, Ds)

    return (ra*sc_kPa*(Rn - QG) + roa*cp*Ds_kPa - ra*Q_LE*(sc_kPa + gamma) ) / ( gamma * Q_LE )

def compute_missing_variables(Ts, Ds):
    Ts_kelvin = Ts + zero_celsius_in_kelvin

    # convert units
    Ds_kPa = Ds / 1000

    #print(f"computing Q_LE(rs={rs}, ra={ra}, Rn={Rn}, QG={QG}, Ds={Ds}, Ts={Ts})")
    lambda_ = b0 * np.power(Ts_kelvin / (Ts_kelvin + b1), 2) # J/kg correct
    gamma = specific_heat_air * air_pressure_at_sea_level_kPa / (specific_gravity_water_vapor * lambda_) # correct

    # [Pa] / (constant [J / kg * K] * [K]) = [kg/m³]
    roa = air_pressure_at_sea_level / (287.05 * Ts_kelvin) # kg / m³ probably correct


    sc = compute_sc(Ts) # Pa / C
    sc_kPa = sc / 1000 # kPa / C
    return sc_kPa, roa, Ds_kPa, gamma

# Taken from https://edis.ifas.ufl.edu/publication/AE459
# Ts is in Celsius °C
# Outputs sc in Pa/°C
def compute_sc(Ts):
    b0, b1, b2 = 0.6108, 17.27, 237.3  # empirical coefficients
    saturated_vapour_pressure = 1e3 * b0 * np.exp(b1 * Ts / (b2 + Ts))  # Pa
    return b1 * b2 * saturated_vapour_pressure / np.power(b2 + Ts, 2) # NOTE: remplace np by torch?


#def Q_LE(gsCO2, Tf, Ts, Pre, Pre0, ra, sc, Rn, QG, roa, cp, es, ea, gamma):
#    ... # NOTE: Other version with other dependancies, probably not useful
