# PM Equation dimensional analysis
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