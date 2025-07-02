%%%%%%%%%%%% VEGETATION COMPOSITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Yakutsk Spasskaya Pad larch
%%% Larix Cojanderi  190 yr old forest
% Canopy height 20m
% LAI 1.5-3
% continuous permafrost zone with sandy loam soils
% GPP 635  NPP 489 AGBi 26 gC/m2 yr-1,
% Tower height 32m
% The forest floor is fully covered by dense cowberry (Vaccinium
% vitis-idaea)
%% Root depth 20cm
%leaf fall beginning fluctuated from DOY 240 in 1998 to DOY 265 and it last
%10-20 days


%%% SOIL AND HYDROLOGICAL PARAMETER
cur_dir=cd;
cd(Directory)
%%% Rainfall Disaggregation
%%%%%
a_dis = NaN ;
pow_dis = NaN;
%%%%%%%%%%%%%%%%%%%%
fpr=1;
SvF=1; %% Sky View Factor
SN=0; %% Stream Identifier
Slo_top=0;  %% [fraction dy/dx]
Slo_pot=zeros(1,ms); %% [fraction dy/dx]
Asur = 1./cos(atan(Slo_top)); %% Real Area/Projected Area [m^2/m^2]
Ared = 1; %% 1-Frock, Reduction in Area due to rock content
aR =1; %%% anisotropy ratio
%Kh=Ks*aR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsize=1; %%[m^2];
aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
Kbot = NaN;%0.5; %% [mm/h] Conductivity at the bedrock layer
Krock = NaN; %% [mm/h] Conductivity of Fractured Rock
zatm = 32; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%%% LAND COVER PARTITION
Cwat = 0; Curb = 0 ; Crock = 0;
Cbare = 0; Ccrown = [1.0];
%%%%%%%%%% SOIL INPUT -
Pcla= 0.10;
Psan= 0.65;
Porg= 0.05;
Color_Class = 0;
%%%%%%%%%%%%%%%%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%%%%%%%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s =lan_s*ones(1,ms);
cv_s = cv_s*ones(1,ms);
%%%
SPAR=2; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%%%
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%% Measured %%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%
%%%%%%%%%%%%%
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
clear Oss_Lp Owp_Lp
%%%%%%%%%%%%%
Zs= [0 10 20 50 100 150 200 300 400 500 600 700 800 1000 1250 1500 1750 2000 ]; %% ms+1
Zdes = 10;
Zinf=  10;
Zbio = 250;
if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end
[EvL_Zs]=Evaporation_layers(Zs,Zdes); %%% Evaporation Layer fraction
[Inf_Zs]=Evaporation_layers(Zs,Zinf); %%% Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio); %%% Infiltration Depth Layer fraction
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
In_max_urb=5;
In_max_rock=0.1; %% [mm]
%%%%%%%%%%%%% SNOW PARAMETER
TminS=-0.8;%% Threshold temperature snow
TmaxS= 2.8;%% Threshold temperature snow
ros_max1=550; %600; %%% [kg/m^3]
ros_max2=260; %450; %%% [kg/m^3]
Th_Pr_sno = 10.0; %%% [mm/day] Threshold Intensity of snow to consider a New SnowFall
%%%% ICE Parameter
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
Aice = 0.35; %% [-] Ice albedo
WatFreez_Th = -8; %% [°C] Threshold for freezing water
dz_ice = 0.54; %% [mm/h] Freezing Layer depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExEM = 1.0;
%%%% PARAMETERS VEGETATION
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile
ZR95_H = [300]; %% [mm]
ZR95_L = [0]; %% [mm]
ZR50_H = [NaN];
ZR50_L = [NaN];
ZRmax_H = [NaN];
ZRmax_L = [NaN];
%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
%5 Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_L_In= [0.2]; %%[mm/LAI]
Sp_LAI_H_In= [0.1]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [0.8]; %%[cm]
d_leaf_L= [0.8];  %% [cm]
%%%%%%%% Biochemical parameter
KnitH=[ 0.20]; %%% Canopy Nitrogen Decay
%KnitH=0.5; %%% Canopy Nitrogen Decay
KnitL=[ 0.0];
%%%%%
mSl_H = [0.0]; %% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = [0.0];  % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%------
%%------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_H=[ 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[ 700]; %%[Pa]
a1_H=[ 6.0];%% 9
go_H=[ 0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[ 0.660];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[94];% [88 88]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[ Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[ 1.5]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_L=[  NaN];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[  NaN]; %%[Pa]
a1_L=[  NaN];  %%% [-] WUE parameter
go_L=[ NaN];% [mol / s m^2] minimum Stomatal Conductance
CT_L=[  NaN]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[  NaN];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[  NaN]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[  NaN]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L=[ NaN]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Matric Potential
Psi_sto_00_H = [ -0.8]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H = [  -2.5];%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H = [ -1.0]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H = [ -3.2] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H =[ 10] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  =[1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = [ 15.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [ 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [ -5.0]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [ 150]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata
%% Hydraulic Parameters
Psi_sto_00_L = [ NaN]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_L = [  NaN] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L = [ NaN]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_L = [  NaN] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [ NaN] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [ NaN];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_L = [  NaN] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [  NaN];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [ NaN]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [ NaN]; %%% [kg / m^3 sapwood MPa]


%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);


%%%% Growth Parameters
PsiG50_H= [  -0.8];  %%[MPa]
PsiG99_H= [  -2.5];  %%[MPa]
gcoef_H = [ 3.5]; % [gC/m2 day]
%%------
PsiG50_L= [  NaN ];
PsiG99_L= [  NaN ];
gcoef_L =[ NaN ];  % [gC/m2 day]


%%%%%%%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(3);
[PFT_opt_L(1)]=Veg_Optical_Parameter(0);
OM_H=[1];
OM_L=[1];
%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.020]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [ 26] ; % [38]; %[kgC/kgN ] Leaf Nitrogen Concentration
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
r_H = [ 0.055];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ] -- Bonan 2003
gR_H= [ 0.25]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
aSE_H= [ 1]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_H= [ 1/200];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  [ 36/365] ; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H = [ 0]; % %% [°C] Cold Leaf Shed
drn_H=  [ 1/1100]; %% turnover root  [1/d]
dsn_H= [ 1/950]; % normal transfer rate sapwood [1/d]
age_cr_H= [ 180]; %% [day] Critical Leaf Age
Bfac_lo_H= [ 0.99]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [ -2]; %% Mean Temperature for Leaf onset
Tls_H = [NaN]; %% Mean Temperature for Leaf Shed
PAR_th_H= [ NaN];
dmg_H= [ 30]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [ 0.01];
Trr_H = [ 3.0]; %% Translocation rate [gC /m^2 d]
mjDay_H = [ 250]; %% Maximum Julian day for leaf onset
LDay_min_H =[ 14.8]; %% Minimum Day duration for leaf onset
LtR_H = [ 0.8]; %%% Leaf to Root ratio maximum
Mf_H= [ 1/50]; %% fruit maturation turnover [1/d]
Wm_H= [0] ; % wood turnover coefficient [1/d]
eps_ac_H = [1.0 ]; % %% Allocation to reserve parameter [0-1]
LDay_cr_H = [ 14.8]; %%%  Threshold for senescence day light [
Klf_H =[ 1/10]; %% Dead Leaves fall turnover [1/d]
fab_H = [ 0.80]; %% fraction above-ground sapwood and reserve
fbe_H = [ 0.20]; %% fraction below-ground sapwood and reserve
ff_r_H= [ 0.1]; %% Reference allocation to Fruit and reproduction
[ParEx_H(1)]=Exudation_Parameter(0);
[Mpar_H(1)]=Vegetation_Management_Parameter;
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L = [NaN]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L= [NaN ]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
r_L = [NaN ];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_L= [NaN ]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_L= [2.0 10.0];
aSE_L= [NaN ]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_L= [NaN ];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L =  [NaN ]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_L = [NaN ]; %% [°C] Cold Leaf SLed
drn_L=  [NaN ]; %% turnover root  [1/d]
dsn_L= [NaN ]; % normal transfer rate sapwood [1/d]
age_cr_L= [NaN ]; %% [day] Critical Leaf Age
Bfac_lo_L= [NaN ]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN ]; %% Leaf Shed Water Stress [0-1]
Tlo_L = [NaN ]; %% Mean Temperature for Leaf onset
Tls_L = [NaN ]; %% Mean Temperature for Leaf Shed
PAR_th_L= [NaN ];
dmg_L= [NaN ]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_L = [NaN ];
Trr_L = [NaN ]; %% Translocation rate [gC /m^2 d]
mjDay_L = [NaN ]; %% Maximum Julian day for leaf onset
LDay_min_L =[NaN ]; %% Minimum Day duration for leaf onset
LtR_L = [NaN ]; %%% Leaf to Root ratio maximum
Mf_L= [NaN ]; %% fruit maturation turnover [1/d]
Wm_L= [NaN ] ; % wood turnover coefficient [1/d]
eps_ac_L = [NaN]; %% Allocation to reserve parameter [0-1]
Klf_L =[NaN]; %% Dead Leaves fall turnover [1/d]
fab_L = [NaN ]; %% fraction above-ground sapwood and reserve
fbe_L =[NaN ]; %% fraction below-ground sapwood and reserve
LDay_cr_L = [NaN]; %%%  Threshold for senescence day light [h]
ff_r_L= [NaN];
[ParEx_L(1)]=Exudation_Parameter(0);
[Mpar_L(1)]=Vegetation_Management_Parameter;
%%%%%%%%%%%% PRODUCTIVITY
%%% Maximum Rubisco Capacity
%[Vmax_H]=Maximum_Rubisco_Capacity(Sl_H,PLNR_H,Nl_H);
%[Vmax_L]=Maximum_Rubisco_Capacity(Sl_L,PLNR_L,Nl_L);
Vmax_H = [38];%
Vmax_L = [0];
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%% Initial Condtion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%
LAI_H(1,:)=[ 0.0]; %
B_H(1,:,:)= [ 0 109 114 79 1.1 0 0.4 0 ]; %%
Rrootl_H(1,:)= [ 3000] ;
PHE_S_H(1,:)=[1];
dflo_H(1,:)=[ 0];
AgeL_H(1,:)=[ 0];
e_rel_H(1,:)=[ 1];
hc_H(1,:) =[ 20.0]; %%
SAI_H(1,:) = [ 0.1]; %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LAI_L(1,:)=[0.0 ]; %
B_L(1,:,:)= [0 0 0 0 0 0 0 0]; %% 95 58 120 18
Rrootl_L(1,:)= [ 0] ;
PHE_S_L(1,:)=[ 0];
dflo_L(1,:)=[0 ];
AgeL_L(1,:)=[0 ];
e_rel_L(1,:)=[1 ];
hc_L(1,:) =[0.0]; %% 0.7
SAI_L(1,:) = [0.0 ]; %% 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [1000];
Preserve_H(1,:)= [100];
Kreserve_H(1,:)= [1000];
FNC_H(1,:)=[1];
NupI_H(1,:,:)= [0 0 0 ];
Nreserve_L(1,:)= [1000];
Preserve_L(1,:)= [100];
Kreserve_L(1,:)= [1000];
FNC_L(1,:)= [1];
NupI_L(1,:,:)= [0 0 0 ];
RexmyI(1,:)= [0 0 0];

if OPT_SoilBiogeochemistry == 1

end
%%%
TBio_L=[ 0 ];  %%[ton DM / ha ]
TBio_H=[ 200 ]; %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=[10];  %% [mm/ m2 PFT];
Vl_H=[10];  %% [mm/ m2 PFT];
Vx_L=[10];   %% [mm/ m2 PFT];
Vl_L=[10];   %% [mm/ m2 PFT];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+2;
Tdamp(1)=0;
Tdp(1,:)=  [ -7.1484   -7.0916   -6.9806   -6.7675   -6.4952   -6.2181   -5.8155   -5.3142   -4.8281   -4.3585   -3.9030   -3.4724   -2.8333   -1.8047  -0.6523   -0.0227   -0.0317];
TdpI_H(1,:)= -4.9235;
TdpI_L(1,:)=0;
%%% Snow_alb = soil_alb initial
snow_alb.dir_vis = 0.7;
snow_alb.dif_vis = 0.7;
snow_alb.dir_nir = 0.7;
snow_alb.dif_nir = 0.7;
In_L(1,:)=0; In_H(1,:)=0;
In_urb(1)=0; In_rock(1)= 0;
In_Litter(1)=0;
SP_wc(1)=0 ; %%[mm]
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
O(1,:)= [ 0.0736    0.0738    0.0740    0.0744    0.0750    0.0756    0.0766    0.0779    0.0794    0.0810    0.0827    0.0846    0.0880    0.0960 0.1171    0.2658    0.2097];
Oice(1,:)= [ 0.1431    0.0990    0.1604    0.1044    0.0504    0.0271    0.0678    0.1185    0.1397    0.2031    0.2197    0.3080    0.2512    0.1940 0.1120    0.2503    0.3084];
%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
Vice(1,:) = (Oice(1,:)).*dz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(cur_dir)