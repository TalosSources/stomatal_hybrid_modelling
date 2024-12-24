"""
This file contains implementations of the differentiable pipelines
relating rs and Q_LE.
To include: 
* the normal pipeline, with differentiable PM Equation
* variations of the normal pipeline, including the function mapping the output of the FCN
* the inv pipeline, a FCN computes rs, trained to match the inverted PM rs_hat reconstruction from Q_LE
"""