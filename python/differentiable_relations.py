"""
In this file are differentiable functions corresponding to 
known relationships and formulas that can be used in a 
differentiable pipeline 
"""
import numpy as np


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
Implements the PM-Equation:
    Computes Q_LE as a function of the stomatal resistance rs, and other predictors.
Parameters:
* rs - stomatal resistance
* ra - aerodynamic resistance
* Rn - Net Radiation
* QG - Soil Heat Flux [W/m²]
* Ds - Vapour-Pressure deficit
* Ts - temperature [°C]
"""
def Q_LE(rs, ra, Rn, QG, Ds, Ts):

    sc_kPa, roa, Ds_kPa, gamma = compute_missing_variables(Ts, Ds)
    
    q_le = (sc_kPa * (Rn - QG) + roa*cp*Ds_kPa / ra ) / (sc_kPa + gamma*(1 + rs/ra)) # flux, W/m²

    return q_le

def compute_rs(Q_LE, ra, Ts, Ds, Rn, QG):
    sc_kPa, roa, Ds_kPa, gamma = compute_missing_variables(Ts, Ds)

    return (ra*sc_kPa*(Rn - QG) + roa*cp*Ds_kPa - ra*Q_LE*(sc_kPa + gamma) ) / ( gamma * Q_LE )

def compute_missing_variables(Ts, Ds):
    Ts_kelvin = Ts + zero_celsius_in_kelvin

    # convert units
    Ds_kPa = Ds / 1000

    lambda_ = b0 * np.power(Ts_kelvin / (Ts_kelvin + b1), 2)
    gamma = specific_heat_air * air_pressure_at_sea_level_kPa / (specific_gravity_water_vapor * lambda_)

    # [Pa] / (constant [J / kg * K] * [K]) = [kg/m³]
    roa = air_pressure_at_sea_level / (287.05 * Ts_kelvin) # kg / m³


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

