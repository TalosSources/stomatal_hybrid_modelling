##############################################################
#   Subfunction  Soil-Plant-Atmosphere-Continuum           ###
##############################################################

# imports
from math import *
import numpy as np
import torch

"""
AUTOMATICALLY DIFFERENTIABLE TRANSLATION OF THE MODULE OF THE SAME NAME
IN T&C, FROM MATLAB TO PYTHON
Args:
    * All the args from the matlab module, torch tensors (scalars or batches)
    * rs_model: parametrized model for rs 
        taking as input 7 numbers (An,Pre,Cc,GAM,Ds,Do, Ts)
    * gsCO2_model: parametrized model for gsCO2 
        taking as input 6 numbers (An,Pre,Cc,GAM,Ds,Do)
    If rs_model and gsCO2_model are both None, we use the empirical formulation for rs.
    * Vmax_model: For now, always expected to be None
Output:
    rs : prediction of the stomatal resistance
"""
def photosynthesis_biochemical(Cc,IPAR,Csl,ra,rb,Ts,Pre,Ds,Psi_L,Psi_sto_50,Psi_sto_00,
    CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,gmes,rjv,
    gsCO2_model = None, rs_model=None, Vmax_model = None):

    # CT is a flag (can be a tensor), must be 3 or 4, and the same for all elements
    CT[torch.isnan(CT)] = 3
    assert (CT == 3).all() or (CT == 4).all()
    if CT.dim() > 0:
        CT = CT[0]
    
    Ta=Ts 
    Pre0 = 101325 ## [Pa] 
    Tf = 273.15 ## [K] 

    Pre = Pre*100 ## [Pa]
    IPAR = IPAR*4.57  ### [umolPhotons/s m^2] ### [Dey 2004]

    ra = ra*(0.0224*(Ta+273.15)*Pre0/(Tf*Pre))*1e-6 ## [m^2 s/umolH20]  ### -> 
    rb = rb*(0.0224*(Ta+273.15)*Pre0/(Tf*Pre))*1e-6 ## [m^2 s/umolH20]  ### ->  
    Cc = Cc*1e-6*Pre ## [Pa] - Partial Pressure [Pa*molCO2/molAIR]
    Oa = Oa*1e-6*Pre ## [Pa]
    Csl = Csl*1e-6*Pre ## [Pa] -- Leaf surface CO2 concentration 

    rmes = 1/(1e+6*gmes)  ## [ s m^2 /umolCO2 ] Mesophyl Conductance
    go = go*1e6 ###  [umolCO2 / s m^2] 

    Ts_k = Ts + 273.15 ##[K]
    Tref = 25 + 273.15 ## [K] Reference Temperature
    R =   0.008314 ##  [kJ�/ K mol] Gas Constant


    ANS_TEMP=1 
    if  ANS_TEMP == 1: ## Kattge and Knorr 2007
        Hd =  200 # [kJ/mol]  Deactivation Energy -- constant
        kT =torch.exp(Ha*(Ts_k-Tref)/(Tref*R*Ts_k))*(1+np.exp((Tref*DS - Hd)/(Tref*R)))/(1+torch.exp((Ts_k*DS-Hd)/(Ts_k*R)))
        Vm=Vmax*kT ### ## [umolCO2/ s m^2 ]
        Hd =  200#  [kJ/mol]  Deactivation Energy -- constant
        Ha = 50# [kJ/mol] Activation Energy
        DS = 0.646# [kJ / mol K]  entropy factor
        kT =torch.exp(Ha*(Ts_k-Tref)/(Tref*R*Ts_k))*(1+np.exp((Tref*DS - Hd)/(Tref*R)))/(1+torch.exp((Ts_k*DS-Hd)/(Ts_k*R)))
        Jmax = Vmax*rjv  ## [umol electrons/ s m^2 ]
        Jm= Jmax*kT #### [umol electrons/ s m^2 ]
    else: ### Bernacchi et al., 2001 2003  Bonan et al., 2011
        Hd = 149 ##  [kJ/mol]  Deactivation Energy -- constant
        Ha = 65.33 ##  [kJ/mol] Activation Energy - Plant Dependent
        DS = 0.485 ## [kJ / mol K]  entropy factor - Plant Dependent
        kT =torch.exp(Ha*(Ts_k-Tref)/(Tref*R*Ts_k))*(1+np.exp((Tref*DS - Hd)/(Tref*R)))/(1+torch.exp((Ts_k*DS-Hd)/(Ts_k*R)))
        Vm=Vmax*kT ### ## [umolCO2/ s m^2 ]
        Hd = 152 ## [kJ/mol]  Deactivation Energy -- constant
        Ha = 43.5 ## [kJ/mol] Activation Energy
        DS = 0.495 ## [kJ / mol K]  entropy factor
        kT =torch.exp(Ha*(Ts_k-Tref)/(Tref*R*Ts_k))*(1+np.exp((Tref*DS - Hd)/(Tref*R)))/(1+torch.exp((Ts_k*DS-Hd)/(Ts_k*R)))
        Jmax = Vmax*rjv  ## [umol electrons/ s m^2 ]
        Jm= Jmax*kT #### [umol electrons/ s m^2 ]

    Ha = 53.1 ## [kJ/mol] Activation Energy
    DS = 0.490 ##  [kJ / mol K]  entropy factor
    Hd = 150.65 ## [kJ/mol]  Deactivation Energy

    TPU25=0.1182*Vmax ## [umolCO2/ s m^2 ]
    kT =torch.exp(Ha*(Ts_k-Tref)/(Tref*R*Ts_k))*(1+np.exp((Tref*DS - Hd)/(Tref*R)))/(1+torch.exp((Ts_k*DS-Hd)/(Ts_k*R)))
    TPU=TPU25*kT ##  [umolCO2/ s m^2 ]

    if CT==4: 
        s1=0.3 ## [1/K] 
        s3=0.2 # [1/K] ## 0.3 Cox2001
        Tup = 40 #[�C]
        Tlow = 15 ##[�C]
        f1T= 1/(1 +torch.exp(s1*(Ts - Tup))) ### Temperaure Function 1 for Maximum Rubisco Capacity
        f2T= 1/(1 +torch.exp(s3*(Tlow - Ts)))### Temperaure Function 2 for Maximum Rubisco Capacity
        fT = 2.0**(0.1*(Ts-25)) 
        Vm = Vmax*fT*f1T*f2T ## [umolCO2/ s m^2 ]
        ke25 = 20000*Vmax 
        ke = ke25*fT 

    ANSG = 2 
    if ANSG == 0: 
        fT = 0.57**(0.1*(Ts-25))
        GAM = Oa/(2*2600*fT)*(CT==3) ## [Pa]
    if ANSG == 1: 
        G0= 34.6### 28  [umolCO2/molAIR]
        G1 = 0.0451 # 0.0509
        G2= 0.000347 # 0.001
        T0 = 293.2 ##[K]
        GAM = G0*(1 +G1*(Ts+273.15 -T0) + G2*(Ts+273.15-T0)**2) ### [umolCO2/molAIR] -- CO2 Compensation point - Leuning 1995
        GAM = GAM*1e-6*Pre ## [Pa] - Partial Pressure [Pa*molCO2/molAIR]
    if ANSG == 2: 
        Ha = 37.83   ## [kJ/mol] Activation Energy 
        kT= torch.exp(Ha*(Ts_k-Tref)/(Tref*R*Ts_k))
        GAM25 = 42.75 ## ##[umol / mol]
        GAM25 = GAM25*1e-6*Pre ##[Pa]
        GAM = GAM25*kT ## [Pa] Michaelis-Menten Constant for C0_2


    if CT == 3:
        Ha =79.43  ## [kJ/mol] Activation Energy
        Kc25= 404.9 ##[umol / mol]
        Kc25= Kc25*1e-6*Pre ##[Pa]
        kT= torch.exp(Ha*(Ts_k-Tref)/(Tref*R*Ts_k))
        Kc = Kc25*kT ## [Pa] Michaelis-Menten Constant for C0_2
        ###
        Ha = 36.38  ## [kJ/mol] Activation Energy
        Ko25 = 278.4  ##[mmol / mol]
        Ko25 = Ko25*1e-3*Pre ##[Pa]
        kT= torch.exp(Ha*(Ts_k-Tref)/(Tref*R*Ts_k))
        Ko = Ko25*kT ## [Pa] Michaelis-Menten Constant for O_2

    if CT == 3:
        Ha = 46.39
        DS = 0.490 
        Hd = 150.65 
        Rdark25 = 0.015*Vmax 
        kT =torch.exp(Ha*(Ts_k-Tref)/(Tref*R*Ts_k))*(1+np.exp((Tref*DS - Hd)/(Tref*R)))/(1+torch.exp((Ts_k*DS-Hd)/(Ts_k*R)))
        Rdark = Rdark25*kT
    elif CT ==4:
        fT = 2.0**(0.1*(Ts-25))
        fT3= 1/(1 +torch.exp(1.3*(Ts-55)))### Temperaure Function 3 for Respiration
        Rdark25 = 0.025*Vmax 
        Rdark = Rdark25*fT*fT3 ###  ### [umolCO2/ s m^2 ] ## Leaf Maintenace Respiration / Dark Respiration


    Q = FI*IPAR ##  [umolCO2/ s m^2 ] ## Light Absorbed by Photosystem II in CO2 units 
    d1= 0.7 
    d2= -(Q + Jm/4) 
    d3= Q*Jm/4   ## d1 = 0.95 Leuning 1995 d1 = 0.7 Bonan et al., 2011  
    J=torch.min((-d2+torch.sqrt(d2**2-4*d1*d3))/(2*d1),(-d2-torch.sqrt(d2**2-4*d1*d3))/(2*d1)) # 
    if CT == 3:
        JC = Vm*(Cc -GAM)/(Cc + Kc*(1+Oa/Ko)) ### Gross Assimilation Rate Limited by Rubisco# [umolCO2/ s m^2 ]
        JL = (J)*(Cc -GAM)/(Cc + 2*GAM) ### Gross Assimilation Rate Limited by Light # [umolCO2/ s m^2 ]
        JE = 3*TPU ## Gross Assimilation Rate Limited by Export # [umolCO2/ s m^2 ]
    elif CT==4:
        JC = Vm ### Gross Assimilation Rate Limited by Rubisco# [umolCO2/ s m^2 ]
        JL = Q ### Gross Assimilation Rate Limited by Light # [umolCO2/ s m^2 ]
        JE = ke*Cc/Pre 

    if CT==3:
        b1= 0.98 
        b2= -(JC+JL) 
        b3= JC*JL   ## b1 =0.98  Collatz 1991  Selleers1996b Bonan et al 2011   b1 = 0.83 Cox 1998
    elif CT==4:
        b1= 0.80 
        b2= -(JC+JL) 
        b3= JC*JL   ##  b1 = 0.8 Bonan et al 2011  for C4  0.83 Cox 2001

    JP=torch.min((-b2+torch.sqrt(b2**2-4*b1*b3))/(2*b1),(-b2-torch.sqrt(b2**2-4*b1*b3))/(2*b1))
    

    if CT == 3:
        c1 = 0.95 
        c2 = -(JP +JE) 
        c3= JP*JE  ### c1 =0.95  Collatz 1991 Sellers1996b Bonan et al 2011   c1 =0.90 Cox1998
    elif CT == 4:
        c1 = 0.95 
        c2 = -(JP +JE) 
        c3= JP*JE  ### c1=0.95 Bonan et al 2011 for C4  0.93 Cox 2001

    A=torch.min((-c2+torch.sqrt(c2**2-4*c1*c3))/(2*c1),(-c2-torch.sqrt(c2**2-4*c1*c3))/(2*c1))

    Rgsws=0.02
    p2= log((1 -Rgsws)/Rgsws)/(Psi_sto_00 - Psi_sto_50)## [1/MPa]
    q2=-p2*Psi_sto_50 ##[-]
    Rgsw = 1/(1+torch.exp(p2*Psi_L+q2)) ## [fraction]
    fO = (1-Rgsw) ######
    fO[fO>1]=1  
    fO[fO<0]=0 

    if CT == 3:
        Jfe = A*(Cc + 2*GAM)/(Cc -GAM)
    elif CT==4:
        Jfe= A ## [umolCO2/ s m^2 ]

    fiP0= FI*4 ### [umol Electrons/ umolPhotons]
    fiP = fiP0*Jfe/Q  ## [0.4 max - stress decrease ]  

    dls=1-fiP/fiP0 ## degree of light saturation 

    kf =0.05 
    kd = torch.max(0.03*Ts+0.0773,torch.tensor(0.087)) 
    kn = (6.2473*dls-0.5944)*dls 

    fiF = kf/(kf+kd+kn)*(1-fiP)  # [umol Electrons/ umolPhotons]
    SIF = IPAR*fiF # ### [umol electrons/s m^2]

    k=0.0375*Vmax +8.25  ## [umol m-2 s-1 / W m-2 sr-1 um-1]  
    F755nm =SIF/k ## [W m-2 sr-1 um-1] 

    A = A*fO ## Gross Assimilation Rate [umolCO2/ s m^2 ]

    An = A - Rdark # ## Net Assimilation Rate # [umolCO2/ s m^2 ]

    """
    2 Options: 
        1. simply predict rs directly, don't predict gsCO2 first.
        we rescale the predictors and output by their expected magnitude, for more homogeneous features.
        2. predict gsCO2, and then map it to rs using physical laws.
    """
    if rs_model is not None:
        predictors = torch.stack([An * 1e2,Pre * 1e-4,Cc * 1e-1,GAM * 1e1,Ds * 1e-1,Do * 1e-2, Ts], dim=-1)
        model_output = rs_model(predictors).squeeze()

        # tests to make the model_output behave better with gradient descent

        rs_small = torch.nn.functional.softplus(model_output) # 0 at -inf, ~x fo~ large x, analytic 
        rs = rs_small * 1e2

        return rs
    else:
        if gsCO2_model is not None:
            predictors = torch.stack([An * 1e2,Pre * 1e-4,Cc * 1e-1,GAM * 1e1,Ds * 1e-1,Do * 1e-2], dim=-1)
            model_output = gsCO2_model(predictors).squeeze() # We are unfortunately always truncated. NOTE: Does this impact the gradient?
            model_output = torch.nn.functional.softplus(model_output)
            gsCO2 = go + a1 * model_output
        else:
            # EMPIRICAL MODEL
            gsCO2 = go + a1*An*Pre/((Cc-GAM)*(1+Ds/Do)) ###  [umolCO2 / s m^2] -- Stomatal Conductance
        
        # For differentiability, we could use something like softplus centered on go
        gsCO2 = torch.max(gsCO2, go)

        # Convert the stomatal conductance to the stomatal resistance
        rsCO2=1/gsCO2 ### [ s m^2 / umolCO2 ] Stomatal resistence or Canopy 
        rsH20 = (rsCO2/1.64)*(1e6) ### [ s m^2 / molH20 ] Stomatal resistence or canopy 
        rs = rsH20*(Tf*Pre)/(0.0224*(Ts+273.15)*Pre0) ## [s/m]  Stomatal resistence or Canopy [convert stomatal resistence in terms of water volumes into another unit]
        return rs
