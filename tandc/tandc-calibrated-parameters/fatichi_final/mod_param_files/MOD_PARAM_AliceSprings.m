%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOIL AND HYDROLOGICAL PARAMETER %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mulga (Acacia) woodland 
%% (Acacia aneura) canopy, which is 6.5m tall on average. 
%%% Evergreen
%discontinuous canopy of short (3–7 m), evergreen Acacia trees (Acacia aptaneura and Acacia aneura)
%with an understorey of shrubs, herbs and grasses (mix of C3 and C4, shurbs e.g., Spinifex grass Triodia sp C4) that are conditionally active depending 
%upon moisture avail-ability and season
%Acacia is 74.5% of the land area in the Mulga woodland  basal area is 8 m2ha the understory  dominated by C4 warm-season grasses
%%%% GPP= 320 gC/m2  LAI = 0.3 but 0.8 in 2011 Ta = 22.8 C  Pr 331 mm/yr  ET =263

%Bowenratio (H/LE) was  37.5 (range 0.78–408) in the Mulga woodland 

%%% sandy loam (74:11:15  sand:silt:clay) Soil organic matter is lessthan 1%  bulk density 1.69 g/cm3
%%% sandy clay overlying a 49m deep water table
%% Tower height  windspeed and wind direction (9.25m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
cur_dir=cd;
cd(Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Basic Grid Cell Information - Not relevant for plot scale 
fpr=1;
SvF=1; %% Sky View Factor
SN=0; %% Stream Identifier
Slo_top=0;  %% Slope [fraction dy/dx]
Slo_pot=zeros(1,ms); %% Slope of hydraulic head [fraction dy/dx]
Asur = 1./cos(atan(Slo_top)); %% Real Area/Projected Area [m^2/m^2]
Ared = 1; %%% Reduction due to soil rock content 
aR =1; %%% anisotropy ratio %Kh=Ks*aR;
cellsize=1; %%[m^2];
aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
%%%%%%%%
Kbot = NaN; %% [mm/h] Conductivity at the bedrock layer
Krock = NaN; %% [mm/h] Conductivity of Fractured Rock
zatm =9.25; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%
%%%% LAND COVER and VEGETATION COMPOSITION 
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare = 0.3; Ccrown = [0.18 0.29 0.23];
%%% Rainfall disaggrgation information 
a_dis = NaN ;
pow_dis = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOIL PARAMETERS 
Pcla= 0.15;
Psan= 0.74;
Porg= 0.001;
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
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%%%%%%%%
Osat=Osat*ones(1,ms);
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%
%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
clear Oss_Lp Owp_Lp
Oice = 0;
%%%%%%%%%%%%%
Zs= [0    10    20    50   100  150  200   300   400   600   800  1000  1200  1500  1750  2000  2250  2500  3000 3500 4000 4500 5000 5500 6000]; %% ms+1
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   SNOW PARAMETER
TminS=-1.1;%% Threshold temperature snow
TmaxS=2.5;%% Threshold temperature snow
ros_max1=520; %%% [kg/m^3]
ros_max2=320; %%% [kg/m^3]
Th_Pr_sno = 3; %%% [mm/h] Threshold Intensity of snow to consider a New SnowFall
%%%% ICE Parameter
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
Aice = 0.28; %% [-] Ice albedo
WatFreez_Th = -8; %% [°C] Threshold for freezing water
dz_ice = 0.45; %% [mm/h] Freezing Layer depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%  ROOT PARAMETER 
ExEM = 0; %% Fraction of ectomychorrizae per area 
CASE_ROOT= 3;  %%% Type of Root Profile
ZR95_H = [3000    0    0]; %% [mm]
ZR95_L = [0  400 400]; %% [mm]
ZR50_H = [800  0  0];
ZR50_L = [0 200  200];
ZRmax_H = [NaN  NaN  NaN];
ZRmax_L = [NaN  NaN  NaN];
%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);
%% INTERCEPTION PARAMETERS 
In_max_urb= 5; %% [mm]
In_max_rock= 0.1; %% [mm]
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
%%% Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sp_SN_In= 5.9; %% [mm/LAI]
%%%%%%%%% Interception Parameter
Sp_LAI_H_In= [0.1    0.2    0.2]; %%[mm/LAI]
Sp_LAI_L_In= [0.2   0.2     0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [0.3     2      2]; %%[cm]
d_leaf_L= [2       2     0.8];  %% [cm]
%% Veg Biochemical parameter
KnitH=[0.25     NaN       NaN]; %%% Canopy Nitrogen Decay
KnitL=[NaN     0.4       0.15];
%%%%%
mSl_H = [0  NaN  NaN];%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = [NaN    0    0]; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%
%%------
%%------
FI_H=[0.081     NaN    NaN];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[1000   NaN   NaN]; %%[Pa] 
a1_H=[5  NaN  NaN];  %%% 6 [-] WUE parameter 
go_H=[0.01   NaN     NaN];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3  NaN  NaN]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649      NaN     NaN];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[72  NaN  NaN]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf   NaN    NaN]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[1.9      NaN   NaN]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------
FI_L=[NaN    0.04    0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[NaN  1500   1200]; %%[Pa] 
a1_L=[NaN   5   6];  %%% [-] WUE parameter 
go_L=[NaN   0.01    0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_L=[NaN   4   3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[NaN    0.649     0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[NaN  72   72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[NaN  Inf Inf ]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L=[NaN     1.9    2.1]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%
Vmax_H = [38  0   0]; % [umol CO2 /m2 s] - Maximum Rubisco Capacity 
Vmax_L = [0  30  55]; % [umol CO2 /m2 s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hydraulic Parameters
Psi_sto_00_H = [-0.6           NaN           NaN]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H = [-7.5           NaN           NaN] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H = [-2.0           NaN           NaN]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H = [-8  NaN  NaN] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H = [10  NaN  NaN] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200   NaN   NaN];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_H = [8  NaN  NaN] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000    NaN    NaN];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-10  NaN  NaN]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [80  NaN  NaN]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata
Psi_sto_00_L= [NaN          -0.5          -0.8]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_L = [NaN          -3.5          -4.8] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L = [NaN          -0.8          -1.1]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_L = [NaN            -5.0          -6.5] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [NaN    5    8] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [NaN  1200  1200];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_L = [NaN    6    6] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [NaN  80000  80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [NaN           -10          -9.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [NaN  150  150]; %%% [kg / m^3 sapwood MPa]
%% Growth Parameters
PsiG50_H= [-1.0           NaN           NaN];  %%[MPa]
PsiG99_H= [-7.5           NaN           NaN];  %%[MPa]
gcoef_H = [3.5           NaN           NaN]; % [gC/m2 day]
%%------
PsiG50_L= [NaN          -3.5          -4.8];
PsiG99_L= [NaN            -5.0          -6.5];
gcoef_L = [NaN           3.5           3.5]; % [gC/m2 day]
%%%
%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(1);
[PFT_opt_H(2)]=Veg_Optical_Parameter(0); 
[PFT_opt_H(3)]=Veg_Optical_Parameter(0);
[PFT_opt_L(1)]=Veg_Optical_Parameter(0);
[PFT_opt_L(2)]=Veg_Optical_Parameter(14); 
[PFT_opt_L(3)]=Veg_Optical_Parameter(13); 

OM_H=[1  NaN  NaN];
OM_L=[NaN    1    1];
%%%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  VEGETATION PART %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% High Vegetation 
aSE_H=  [0  NaN  NaN]; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_H = [0.008       NaN      NaN]; %  [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_H = [50  NaN  NaN]; %[kgC/kgN ] Leaf Nitrogen Concentration
r_H =  [0.042      NaN      NaN];  %% respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_H=  [0.25     NaN        NaN];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
dd_max_H= 1./[365  NaN  NaN];%%%  [1/d]  death maximum for drought
dc_C_H =  1./[182.5           NaN           NaN]; %%[Factor of increasing mortality for cold]
Tcold_H = [3  NaN  NaN]; %% [°C] Cold Leaf Shed
drn_H=  1./[900  NaN  NaN]; %% turnover root  [1/d]
dsn_H= 1./[1000  NaN  NaN]; % normal transfer rate sapwood [1/d]
age_cr_H= [980  NaN  NaN]; %% [day] Critical Leaf Age
Trr_H = [0.4     NaN     NaN]; %% Translocation rate [gC /m^2 d]
LtR_H = [0.9     NaN     NaN]; %%% Leaf to Root ratio maximum
Mf_H= 1./[80  NaN  NaN]; %% fruit maturation turnover [1/d]
Wm_H= [0  NaN  NaN] ; % wood turnover coefficient [1/d]
eps_ac_H = [0.8        NaN        NaN]; %% Allocation to reserve parameter [0-1]
Klf_H = 1./[50  NaN  NaN]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.75         NaN       NaN]; %% fraction above-ground sapwood and reserve
fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1     NaN      NaN]; %% Reference allocation to Fruit and reproduction
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
[Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
[Stoich_H(3)]=Veg_Stoichiometric_Parameter(Nl_H(3));
%%%% Phenology 
Bfac_lo_H= [0.99       NaN      NaN]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN  NaN  NaN] ;  % Not-used 
Tlo_H = [8  NaN  NaN]; %% Mean Temperature for Leaf onset
Tls_H = [NaN  NaN  NaN]; %%% Not-used 
PAR_th_H= [NaN  NaN  NaN]; %% Light Phenology Threshold 
dmg_H= [15  NaN  NaN]; %%%  Day of Max Growth
LAI_min_H = [0.001     NaN        NaN];
mjDay_H = [-1  NaN  NaN]; %% Maximum Julian day for leaf onset
LDay_min_H =[10.6      NaN        NaN]; %% Minimum Day duration for leaf onset
LDay_cr_H = [10.6       NaN       NaN]; %%%  Threshold for senescence day light [h]
[ParEx_H(1)]=Exudation_Parameter(0);
[Mpar_H(1)]=Vegetation_Management_Parameter;
[ParEx_H(2)]=Exudation_Parameter(0);
[Mpar_H(2)]=Vegetation_Management_Parameter;
[ParEx_H(3)]=Exudation_Parameter(0);
[Mpar_H(3)]=Vegetation_Management_Parameter;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Low Vegetation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
aSE_L=  [NaN    2    2]; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_L = [NaN     0.010   0.012]; %  [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_L = [NaN   42   50]; %[kgC/kgN ] Leaf Nitrogen Concentration
r_L =  [NaN     0.055    0.060];  %% respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_L=  [NaN     0.25      0.25];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
dd_max_L= 1./[NaN   10   18];%%%  [1/d]  death maximum for drought
dc_C_L =  1./[NaN      55      55]; %%[Factor of increasing mortality for cold]
Tcold_L = [NaN    2    0]; %% [°C] Cold Leaf Shed
drn_L=  1./[NaN  1290  1690]; %% turnover root  [1/d]
dsn_L= 1./[NaN  365  365]; % normal transfer rate sapwood [1/d]
age_cr_L= [NaN  180  180]; %% [day] Critical Leaf Age
Trr_L = [NaN       3.0      3.0]; %% Translocation rate [gC /m^2 d]
LtR_L = [NaN        0.4        0.5]; %%% Leaf to Root ratio maximum
Mf_L= 1./[NaN   50   50]; %% fruit maturation turnover [1/d]
Wm_L= [NaN    0    0] ; % wood turnover coefficient [1/d]
eps_ac_L = [NaN      1     1]; %% Allocation to reserve parameter [0-1]
Klf_L = 1./[NaN   20   20]; %% Dead Leaves fall turnover [1/d]
fab_L = [NaN    0    0]; %% fraction above-ground sapwood and reserve
fbe_L = 1-fab_L; %% fraction below-ground sapwood and reserve
ff_r_L= [NaN           0.1           0.1]; %% Reference allocation to Fruit and reproduction
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
[Stoich_L(2)]=Veg_Stoichiometric_Parameter(Nl_L(2));
[Stoich_L(3)]=Veg_Stoichiometric_Parameter(Nl_L(3));
%%%% Phenology 
Bfac_lo_L= [NaN          0.99          0.99]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN  0.5  0.5] ; % Not-used 
Tlo_L = [NaN    8    8]; %% Mean Temperature for Leaf onset
Tls_L = [NaN  NaN  NaN]; %% Not-used 
PAR_th_L= [NaN  NaN  NaN]; %% Light Phenology Threshold 
dmg_L= [NaN   20   20]; %%%  Day of Max Growth
LAI_min_L = [NaN          0.01          0.01];
mjDay_L = [NaN   -1   -1]; %% Maximum Julian day for leaf onset
LDay_min_L =[NaN   10   10]; %% Minimum Day duration for leaf onset
LDay_cr_L = [NaN   10   10]; %%%  Threshold for senescence day light [h]
[ParEx_L(1)]=Exudation_Parameter(0);
[Mpar_L(1)]=Vegetation_Management_Parameter;
[ParEx_L(2)]=Exudation_Parameter(0);
[Mpar_L(2)]=Vegetation_Management_Parameter;
[ParEx_L(3)]=Exudation_Parameter(0);
[Mpar_L(3)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initial Condtion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%%%%%%%%%%%%%%%
%%% B1 Leaves - Grass  %%% B2 Sapwood  %%% B3 Fine Root  %%% B4 Carbohydrate Reserve
%%% B5 Fruit and Flower %%% B6 Heartwood - Dead Sapwood %%% B7 Leaves - Grass -- Standing Dead
%%%%%%%%%%%%%%
LAI_H(1,:)=[0.39  0  0]; %
B_H(1,:,:)= [59 41 70 288 1 0 5 0 ;0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0]; %% 
Rrootl_H(1,:)= [1073  0  0] ;
PHE_S_H(1,:)=[3  0  0];
dflo_H(1,:)=[3649  0  0];
AgeL_H(1,:)=[1168  0  0];
e_rel_H(1,:)=[1  1  1];
hc_H(1,:) =[6.0  0  0]; %%
SAI_H(1,:) = [0.1  0  0]; %% 
%%%%%%%%%%%%%%%%%%%%%
LAI_L(1,:)=[0  0  0]; %
B_L(1,:,:)= [0 0 0 0 0 0 0 0 ; 33 0 83 35 1 0 1 0  ;5 0 24 20 1 0 1 0]; %% 
Rrootl_L(1,:)= [0  578 609] ;
PHE_S_L(1,:)=[0  2  2];
dflo_L(1,:)=[0  15  15];
AgeL_L(1,:)=[0  8  12];
e_rel_L(1,:)=[1  1  1];
hc_L(1,:) =[0  0.2  0.15]; %%
SAI_L(1,:) = [0  0.001  0.001]; %% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%
PARI_H(1,:,:)= [0 0 0 ; 0 0 0 ; 0 0 0]; 
NBLI_H(1,:)=   [0 0 0];
PARI_L(1,:,:)=[0 0 0 ; 0 0 0 ; 0 0 0]; 
NBLI_L(1,:)=[0 0 0];
%%%%%%%%%%
Nreserve_H(1,:)= [1000  1000  1000];
Preserve_H(1,:)= [1000  1000  1000];
Kreserve_H(1,:)= [1000  1000  1000];
FNC_H(1,:)=[1  1  1];
NupI_H(1,:,:)= [0 0 0 ;0 0 0 ;0 0 0];
Nreserve_L(1,:)= [1000  1000  1000];
Preserve_L(1,:)= [1000  1000  1000];
Kreserve_L(1,:)= [1000  1000  1000];
FNC_L(1,:)=[1  1  1];
NupI_L(1,:,:)= [0 0 0 ;0 0 0 ;0 0 0];
RexmyI(1,:)= [0 0 0];
%%%
%%%%
%%%%
if OPT_SoilBiogeochemistry == 1
    %%%%%
end
%%%
TBio_L=[1  1  1];  %%[ton DM / ha ]
TBio_H=[1  1  1]; %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=[1  1  1];  %% [mm/ m2 PFT];
Vl_H=[1  1  1];  %% [mm/ m2 PFT];
Vx_L=[1  1  1];   %% [mm/ m2 PFT];
Vl_L=[1  1  1];   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%%%% Initialization OTHER
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+1;
Tdamp(1)=Ta(1);
Tdp(1,:)= Ta(1)*ones(1,ms);
Sdp(1,:)= (0.001*dz).*cv_s.*(Tdp(1,:)); %% [J °C/m^2 K]
TdpI_H(1,:)=29;
TdpI_L(1,:)=29;
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
O(1,:)=[     0.1932    0.1972    0.2032    0.2112    0.2168    0.2194    0.2182    0.1739    0.0978    0.0664    0.0664    0.107,...
    0.1089    0.0930    0.0985    0.1028    0.1058    0.1088   0.1028    0.1058    0.1088   0.1028    0.1058    0.1088]; 

%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%%%%%%%%%%%%%%%%%
cd(cur_dir)
%%%%%%%%%%%%%%%%%
