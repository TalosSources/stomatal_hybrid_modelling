%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stipa tenacissima  perennial tussock grass, highly adapted to water stress conditions and has long, narrow leaves,
%a large root:shoot ratio, and shallow roots that spread from both seeds
%and stolons. 
%vegetation cover varying between 23.1 

%“Amoladeras” is dominated by bare soil, gravel and rock (43%),
%with a 23% of Stipa tenacissima. All of these sites are characterized
%by summer droughts and asynchronous annual patterns of sunlight and precipitation.

%the “degraded” site  (grazed) 

%Clay 14.6 Sand 58.4  Soil depth 10 cm 
%Pr 220 
%height vegetation between 60 and 100 cm

%%%%%%%%%%%%%%%%%
cur_dir=cd;
cd(Directory)
%%% Rainfall Disaggregation 
a_dis = NaN ; 
pow_dis = NaN; 
%%%%% General Topographic Parameters 
fpr=1; 
SvF=1; %% Sky View Factor  
SN=0; %% Stream Identifier 
Slo_top=0;  %% [fraction dy/dx]
Slo_pot=zeros(1,ms); %% [fraction dy/dx]
Asur = 1./cos(atan(Slo_top)); %% Real Area/Projected Area [m^2/m^2] 
Ared = 1; %% 1-Frock  Reduction for Rock Content 
aR =1; %%% anisotropy ratio 
%Kh=Ks*aR;
cellsize=1; %%[m^2]; 
aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght  
%%%%%% Hydraulic Subsurface 
Kbot = NaN; %% [mm/h] Conductivity at the bedrock layer 
Krock = NaN; %% [mm/h] Conductivity of Fractured Rock
%%%% Reference height Measurements 
zatm = 2.5; %% Reference Height
%%%%%%%% COMPOSITION 
%%%% VEG. SPECIES  --- Low Grasses -- High Decidous      
%%%% LAND COVER PARTITION
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare = 0.77; Ccrown = [0.15 0.08];
%%%%%%%%%% SOIL INPUT %%% 
Pcla= 0.146; 
Psan= 0.584;
Porg= 0.005; 
Color_Class = 0;  
%%%%%%%%%%%%%%%%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%%%%%%%%
rsd=rsd*ones(1,ms); 
lan_dry=lan_dry*ones(1,ms); 
lan_s =lan_s*ones(1,ms); 
cv_s = cv_s*ones(1,ms);
%%%%%%%%%%%%%%%
SPAR=2; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls 
%%%%%%%%%%%%%%%%
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
%%%%%%%%%%%%%%
clear Oss_Lp Owp_Lp 
%%%%%%%%%%%%%%
Oice = 0;
%%%%%%%%%%%%%
Zs = [ 0 10 20 50 75 100 150 200 250]; %%% [ms+1]
Zdes = 10;
Zinf=  10; 
Zbio = 250; 
if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return 
end
[EvL_Zs]=Evaporation_layers(Zs,Zdes); %%% Evaporation Layer fraction 
[Inf_Zs]=Evaporation_layers(Zs,Zinf); %%% Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio);
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
Dz=zeros(1,ms); 
for i = 1:ms
    if i>1
        Dz(i)= (dz(i)+ dz(i-1))/2; %%% Delta Depth Between Middle Layer  [mm]
    else
        Dz(i)=dz(1)/2; %%% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end
%%%%%%%%%%%%%%%%% OTHER PARAMETER
In_max_urb=5; In_max_rock=0.1; %% [mm]
%%%%%%%%%%%%% SNOW PARAMETER
TminS=-0.8;%% Threshold temperature snow
TmaxS= 2.8;%% Threshold temperature snow
ros_max1=580; %600; %%% [kg/m^3]
ros_max2=300; %450; %%% [kg/m^3]
Th_Pr_sno = 8.0; %%% [mm/day] Threshold Intensity of snow to consider a New SnowFall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ICE Parameter 
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice 
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice 
Aice = 0.28; %% [-] Ice albedo 
WatFreez_Th = -8; %% [°C] Threshold for freezing lake water
dz_ice = 0.45; %% [mm / h] Water Freezing Layer progression without snow-layer  
%%%%%%%%%%%%%%%%%%%% 
ExEM = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile 
%%%%%
ZR95_H = [0 0]; %% [mm]
ZR95_L = [150 250]; %% [mm]
ZR50_H =  [NaN NaN];
ZR50_L =  [NaN NaN];
ZRmax_H =  [NaN NaN];
ZRmax_L =  [NaN NaN];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
%5 Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_L_In= [0.2 0.2]; %%[mm/LAI]
Sp_LAI_H_In= [0.2 0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [ 0 0]; %%[cm]
d_leaf_L= [0.7 2.0];  %% [cm]
%%%%%%%% Biochemical parameter
KnitH=  [0.0 0.0]; %%% Canopy Nitrogen Decay
KnitL=  [0.15 0.25]; %%% Canopy Nitrogen Decay
mSl_H = [0.0 0.0];%% [m2 PFT /gC]  Linear increase in Sla with LAI 
mSl_L = [0.0 0.0]; % 0.001; %% [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree 
%%%%%
%%%%% 
%%------   Photosynthesis Parameters 
FI_H=[0.081 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[1000 1000]; %%[Pa]
a1_H=[7 7];
go_H=[0.01 0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3 3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent 
Ha_H =[72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent 
gmes_H=[Inf Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance 
rjv_H=[1.97 1.97]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_L=[0.081 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[1000 1000]; %%[Pa]
a1_L=[7 6]; %% [6 4]
go_L=[0.01 0.01];% % [mol / s m^2] minimum Stomatal Conductance
CT_L=[3 3];  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[0.649 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent 
Ha_L =[72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent 
gmes_L=[Inf Inf];
rjv_L= [1.9 2.0]; % 
%%%%%%%% HYDRAULIC PARAMETERS 
%%% Stomata
Psi_sto_00_H =  [-1.4 -1.5]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H =  [-3.0 -3.5] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H =  [-1.3 -1.5] ;%%[MPa]  Water Potential at 50% loss conductivity
PsiL50_H =  [-4.5 -7.5]; %% [MPa]  Water Potential at 2% loss conductivity
Kleaf_max_H = [ 20 20] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200 1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = [6.0 8.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-6.5 -13.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [80 80]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata %% -1.4 4.0 
Psi_sto_00_L =  [-1.2 -1.2];% [-0.7 -0.7]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_L = [-5.5 -6.0]; %[-2.5 -5.5] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L =  [-1.5 -1.4]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_L =  [-6.0 -7.0] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [ 5 8] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [1200 1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L = [6.0 6.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [-10 -9.0]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [150 150]; %%% [kg / m^3 sapwood MPa]
%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L); 

%%%% Growth Parameters 
PsiG50_H= [-1.4;-0.45];  %%[MPa]
PsiG99_H= [-1.2; -1.2] ; %%[MPa]
gcoef_H = [3.5 3.5]; % [gC/m2 day]
%%------  
PsiG50_L= [-2.5;-1.4];  %%[MPa]
PsiG99_L= [-3.0; -5.5] ; %%[MPa]
gcoef_L = [3.5 3.5]; % [gC/m2 day]
%%%%%%%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(0);
[PFT_opt_H(2)]=Veg_Optical_Parameter(0); 
[PFT_opt_L(1)]=Veg_Optical_Parameter(13); 
[PFT_opt_L(2)]=Veg_Optical_Parameter(9); 
OM_H=[1 1];
OM_L=[1 1];
%%%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.015 0.015]; % 0.018 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [30 30]; %[gC/gN ] Leaf Carbon-Nitrogen ratio 
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
[Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
r_H = [0.030 0.030];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_H= [0.25 0.25]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_H= [2.0 10.0];
aSE_H= [1 1]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops 
dd_max_H= [1/365 1/365]; %%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  [2/365 2/365]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H = [7 7]; %% [°C] Cold Leaf Shed
drn_H=  [1/1095 1/1095]; %% turnover root  [1/d]
dsn_H= [1/365 1/365]; % normal transfer rate sapwood [1/d] 
age_cr_H= [150 150]; %% [day] Critical Leaf Age
Bfac_lo_H= [0.95 0.95]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [12.9 12.9]; %% Mean Temperature for Leaf onset
Tls_H = [NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_H= [-Inf -Inf]; 
dmg_H= [35 35]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.01 0.01];
Trr_H = [3.5 3.5]; %% Translocation rate [gC /m^2 d]
mjDay_H = [180 180]; %% Maximum Julian day for leaf onset
LDay_min_H =[11.0 11]; %% Minimum Day duration for leaf onset 
LtR_H = [1.0 1]; %%% Leaf to Root ratio maximum 
Mf_H= [0 0 ]; %% fruit maturation turnover [1/d]
Wm_H= [0 0] ; % wood turnover coefficient [1/d]
eps_ac_H = [1 1]; %% Allocation to reserve parameter [0-1] 
LDay_cr_H = [12.30 12.3]; %%%  Threshold for senescence day light [h]
Klf_H =[1/15 1/15]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74 0.74]; %% fraction above-ground sapwood and reserve
fbe_H = [0.26 0.26]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1 0.1]; %% Reference allocation to Fruit and reproduction 
[ParEx_H(1)]=Exudation_Parameter(0);
[ParEx_H(2)]=Exudation_Parameter(0);  
[Mpar_H(2)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Grass  Shrub 
Sl_L = [0.012 0.010]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L= [75 40]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
[Stoich_L(2)]=Veg_Stoichiometric_Parameter(Nl_L(2));
%PLNR_L = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_L = [0.030 0.036];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_L= [0.25 0.25]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_L= [2.0 10.0];
aSE_L= [2 0]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_L= [1/12 1/365]; %[1/50 1/50];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L = [7/365 2/365]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_L = [1.0 1.0]; %% [°C] Cold Leaf SLed
drn_L= [1/1600 1/900 ]; %% [1/900 1/1200]; turnover root  [1/d]
dsn_L= [1/365  1/600 ]; % normal transfer rate sapwood [1/d]
age_cr_L= [455 730 ]; %% [day] Critical Leaf Age
Bfac_lo_L= [0.99 0.99]; %% Leaf Onset Water Stress
Bfac_ls_L= [0.15 0.15]; %% Leaf Shed Water Stress [0-1]
Tlo_L = [8.0 11.0]; %% Mean Temperature for Leaf onset 
Tls_L = [NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_L= [-Inf -Inf]; 
dmg_L= [20 15]; %%% 20 -- Tree 30 Grasses Day of Max Growth
LAI_min_L = [0.05 0.001 ];
Trr_L = [1.6 0.4]; %% 1.3 Translocation rate [gC /m^2 d]
mjDay_L = [366 366]; %% Minimum Julian day for leaf onset
LDay_min_L =[9.0 9.0]; %% Day duration for leaf onset
LtR_L =  [0.40 0.5]; %%% Leaf to Root ratio maximum
Mf_L= [1/50 1/50]; %% fruit maturation turnover [1/d]
Wm_L= [0 0] ; % wood turnover coefficient [1/d]
eps_ac_L = [1.0 0.8]; %% [0.5] Allocation to reserve parameter [0-1]
LDay_cr_L = [9.0 9.0] ; %[11.0 11.0]; %%%  Threshold for senescence day light [h] [12.4 12.4]
Klf_L = [1/50 1/50]; % [1/100 1/100];% ; %% Dead Leaves fall turnover [1/d]
fab_L = [0.0 0.75 ]; %% fraction above-ground sapwood and reserve
fbe_L = [1.0 0.25 ]; %% fraction below-ground sapwood and reserve
ff_r_L= [0.1 0.1];
[ParEx_L(1)]=Exudation_Parameter(0); 
[ParEx_L(2)]=Exudation_Parameter(0);  
[Mpar_L(1)]=Vegetation_Management_Parameter; 
[Mpar_L(2)]=Vegetation_Management_Parameter; 
%%%%%%%%%%%%% PRODUCTIVITY
Vmax_H = [0 0]; 
Vmax_L = [60 48]; % Vmax_L = [75 55]; 
%%%%%%%%%%%%%%%%%%%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
 [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end 
Lmax_day = max(L_day); 
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%%%%%%%%%%%%%%%%%%% Initialization VEG 
%%%%%%%%% DORMANT 1 - MAX GROWTH 2 - NORMAL GROWTH 3 - SENESCENCE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%
LAI_H(1,:)=[0.0]; %
B_H(1,:,:)= [0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0]; %%
Rrootl_H(1,:)= [0 0] ; 
PHE_S_H(1,:)=[1 1];
dflo_H(1,:)=[0 0];
AgeL_H(1,:)=[0 0];
e_rel_H(1,:)=[1 1];
hc_H(1,:) =[0 0]; %% 0.7
SAI_H(1,:) = [0 0]; %% 0.1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LAI_L(1,:)=[0.4];
B_L(1,:,:)= [71 0 376 252 4 0 4 0; 
    175 293 351 196 8 0 10 0];
Rrootl_L(1,:) = [2000 2000] ; 
PHE_S_L(1,:)=[3 3];
dflo_L(1,:)=[0 0];
AgeL_L(1,:)=[56 652];
e_rel_L(1,:)=[1 1];
hc_L(1,:) =[0.7 0.4];
SAI_L(1,:) = [0.001 0.1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [0 0 ];
Preserve_H(1,:)= [0 0];
Kreserve_H(1,:)= [0 0];
FNC_H(1,:)=[1 1];
NupI_H(1,:,:)= [0 0 0 ; 0 0 0];
Nreserve_L(1,:)= [1000 1000];
Preserve_L(1,:)= [100 100];
Kreserve_L(1,:)= [100 100];
FNC_L(1,:)=[1 1];
NupI_L(1,:,:)= [0 0 0 ; 0 0 0];
RexmyI(1,:)= [0 0 0];
%%%%%%%%%%%%%%%%%%%%%%
if OPT_SoilBiogeochemistry == 1
end
%%%%
Vx_H=[0 0];  %% [mm/ m2 PFT];
Vl_H=[0 0];  %% [mm/ m2 PFT];
Vx_L=[10 10];   %% [mm/ m2 PFT];
Vl_L=[100 100];   %% [mm/ m2 PFT];
%%%
TBio_L=[1 1];  %%[ton DM / ha ]
TBio_H=[0 0];  %[ton DM / ha ]
%%%%%%%%%%%%%%%%% Initialization OTHER  
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+2;
Tdamp(1)=Ta(1); 
Tdp(1,:)= Ta(1)*ones(1,ms);
TdpI_H(1,:)=15;
TdpI_L(1,:)=15;
Sdp(1,:)= (0.001*dz).*cv_s.*(Tdp(1,:)); %% [J °C/m^2 K]
%%% Snow_alb = soil_alb initial 
snow_alb.dir_vis = 0.2;
snow_alb.dif_vis = 0.2;
snow_alb.dir_nir = 0.2;
snow_alb.dif_nir = 0.2;
In_L(1,:)=0; In_H(1,:)=0;
In_urb(1)=0; In_rock(1)= 0;
In_Litter(1)=0;
SP_wc(1)=0 ; 
In_SWE(1)= 0;
ros(1)= 0;
t_sls(1)= 0;
e_sno(1) = 0.97;
tau_sno(1) = 0;
EK(1)=0; 
WAT(1) = 0;
ICE(1) = 0; 
IP_wc(1)=0; 
ICE_D(1)= 0;
FROCK(1)=0;
Ws_under(1)=1; 
%%%%%%%%%%%%%% Volume [mm]
O(1,:)=  [  0.1786    0.1790    0.1800    0.1817    0.1833    0.1853    0.1873    0.1883];
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
cd(cur_dir)
%%%%%%%%%%%%%%%%%