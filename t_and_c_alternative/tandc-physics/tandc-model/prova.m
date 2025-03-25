load('/home/talos/git_epfl/stomatal_hybrid_modelling/t_and_c_alternative/CO2_Data/Ca_Data.mat')
load('/home/talos/git_epfl/stomatal_hybrid_modelling/t_and_c_alternative/tandc-forcing/Data_AU-ASM_run.mat')
current_directory = cd;
NN= 43824;  %%% time Step 
dt=3600; %%[s] %%% 
dth=1; %%[h]
ms=24; %%% Soil Layer 
cc = 3; %% Crown area 
id_location = 'AliceSprings'; 
NN=length(Date);
N=Latm; 
x1=1; 
x2=NN;  
Date=Date(x1:x2); 
Pr=Pr(x1:x2);
Ta=Ta(x1:x2);
Ws=Ws(x1:x2); ea=ea(x1:x2); Pre=Pre(x1:x2); SAD1=SAD1(x1:x2); 
SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2); 
SAB2=SAB2(x1:x2); N=N(x1:x2);Tdew=Tdew(x1:x2);esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD = PARD(x1:x2); 
t_bef= 0.5; t_aft= 0.5;
Ds=esat-ea; %% [Pa] Vapor Pressure Deficit 
Ds(Ds<0)=0;
d1 = find(abs(Date_CO2-Date(1))<1/36);d2 = find(abs(Date_CO2-Date(end))<1/36);
Ca=Ca(d1:d2); 
clear d1 d2 Date_CO2 
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
Ws(Ws<=0)=0.01; 
[YE,MO,DA,HO,MI,SE] = datevec(Date); 
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO; 
clear YE MO DA HO MI SE
PARAM_IC='/home/talos/git_epfl/stomatal_hybrid_modelling/t_and_c_alternative/tandc-physics/tandc-model/mod_param.m'
PARAM_IC='/home/talos/git_epfl/stomatal_hybrid_modelling/t_and_c_alternative/tandc-physics/tandc-model/mod_param.m'
Directory ='/home/talos/git_epfl/stomatal_hybrid_modelling/t_and_c_alternative/tandc-physics/tandc-model/src'
Directory='/home/talos/git_epfl/stomatal_hybrid_modelling/t_and_c_alternative/tandc-physics/tandc-model/src'
cd(Directory)
MAIN_FRAME ; 
save('/home/talos/git_epfl/stomatal_hybrid_modelling/t_and_c_alternative/tandc-physics-output/fatichi_parameters_ismail/Results_AU-ASM.mat')
cd(current_directory); 
