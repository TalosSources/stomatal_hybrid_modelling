%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The leaf area index (LAI) generally ranges between 0 and 2
% 
% The site is a typical low tree and shrub savanna environment with _3% tree cover (Rasmussen et al., 2011).
%  The most abundant tree species are Balanites aegyptiaca, Acacia tortilis and Acacia Senegal.   Grass species herbaceous species were 
% Zornia latifolia, Aristida adscensionis, Cenchrus biflorus, Dactyloctenium aegyptium, and Eragrostis
% Tremula
% 
% Acacia Senegal - small thorny deciduous tree
% Balanites aegyptiaca – semi-Deciduous 
% 
% average tree height was 5.2 m
% 
% The soil is sandy luvic arenosol with low amounts of organic material and low clay content (clay = 0.35%, silt = 4.61%, and sand = 95.04%). 
% The land is grazed and located within the Centre de Recherches Zootechniques
% 
% (clay = 0.35%, silt = 4.61%, and sand = 95.04%).
% 
% C3/C4 --   Mostly C4 ;;   16% C3 and 84% C4 
% 
% GPP = 880-1300 
% Biomass 206-471 gDM /m2 
% Tower height 9.0 
% 
% 2010-2013 – Pr 650 486 606 355  
% 
% Grazing - grazing is permanent and occurs year-aroud



%%% SOIL AND HYDROLOGICAL PARAMETER %
%%%%%%%%%%%%%%%%%
cur_dir=cd;
cd(Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
zatm =9.0; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%
%%%% LAND COVER and VEGETATION COMPOSITION 
%%%% Acacia Deciduous / C4 grass / C3 Grass 
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare = 0.0;  Ccrown = [0.03 0.82 0.15];
%%% Rainfall disaggrgation information 
a_dis = NaN ;
pow_dis = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOIL PARAMETERS 
Pcla= 0.01;
Psan= 0.95;
Porg= 0.025;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
Zs= [0    10    20    50   100  150  200   300   400   600   800  1000  1200  1500  1750  2000  2250  2500  3000]; %% ms+1
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
CASE_ROOT= 1;  %%% Type of Root Profile
ZR95_H = [600     0     0]; %% [mm]
ZR95_L = [0  500  400]; %% [mm]
ZR50_H = [NaN  NaN  NaN];
ZR50_L = [NaN  NaN  NaN];
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
Sp_LAI_H_In= [0.2         0.2         0.2]; %%[mm/LAI]
Sp_LAI_L_In= [0.2         0.2         0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [2  2  2]; %%[cm]
d_leaf_L= [2           2         0.8];  %% [cm]
%% Veg Biochemical parameter
KnitH=[0.3           NaN           NaN]; %%% Canopy Nitrogen Decay
KnitL=[NaN           0.4          0.3];
%%%%%
mSl_H = [0  NaN  NaN];%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = [NaN    0    0]; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%
%%------
%%------
FI_H=[0.081           NaN           NaN];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[800  NaN  NaN]; %%[Pa] 
a1_H=[7  NaN  NaN];  %%% [-] WUE parameter 
go_H=[0.01           NaN           NaN];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3  NaN  NaN]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.635          NaN           NaN];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[72  NaN  NaN]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf    NaN    NaN]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[2.4  NaN  NaN]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
%%%------
FI_L=[NaN      0.04      0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[NaN  1000   1200]; %%[Pa] 
a1_L=[NaN    4    6];  %%% [-] WUE parameter 
go_L=[NaN          0.01          0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_L=[NaN    4    3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[NaN         0.649         0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[NaN   72   72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[NaN  Inf  Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L=[NaN     1.9      2.1]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%
Vmax_H = [82   0   0]; % [umol CO2 /m2 s] - Maximum Rubisco Capacity 
Vmax_L = [0  40  62]; % [umol CO2 /m2 s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hydraulic Parameters
Psi_sto_00_H = [-1.5           NaN           NaN]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H = [-2.8           NaN           NaN] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H = [-1.8           NaN           NaN]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H = [-2.5           NaN           NaN] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H = [10  NaN  NaN] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200   NaN   NaN];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_H = [5  NaN  NaN] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000    NaN    NaN];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-6  NaN  NaN]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [150  NaN  NaN]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata
Psi_sto_00_L= [NaN          -0.5          -0.8]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_L = [NaN          -2.5          -4.8] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L = [NaN          -1.3          -1.1]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_L = [NaN          -3.5          -3.5] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [NaN    5    8] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [NaN  1200  1200];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_L = [NaN    6    6] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [NaN  80000  80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [NaN           -10          -9.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [NaN  150  150]; %%% [kg / m^3 sapwood MPa]
%% Growth Parameters
PsiG50_H= [-1.5           NaN           NaN];  %%[MPa]
PsiG99_H= [-2.8           NaN           NaN];  %%[MPa]
gcoef_H = [3.5           NaN           NaN]; % [gC/m2 day]
%%------
PsiG50_L= [NaN          -2.5          -3.5];
PsiG99_L= [NaN          -3.5          -4.8];
gcoef_L = [NaN           3.5           3.5]; % [gC/m2 day]
%%%
%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(4);
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
Sl_H = [0.016           NaN           NaN]; %  [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_H = [40  NaN  NaN]; %[kgC/kgN ] Leaf Nitrogen Concentration
r_H =  [0.04           NaN           NaN];  %% respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_H=  [0.25           NaN           NaN];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
dd_max_H= 1./[30  NaN  NaN];%%%  [1/d]  death maximum for drought
dc_C_H =  1./[365  NaN  NaN]; %%[Factor of increasing mortality for cold]
Tcold_H = [5  NaN  NaN]; %% [°C] Cold Leaf Shed
drn_H=  1./[1610  NaN  NaN]; %% turnover root  [1/d]
dsn_H= 1./[1250  NaN  NaN]; % normal transfer rate sapwood [1/d]
age_cr_H= [250  NaN  NaN]; %% [day] Critical Leaf Age
Trr_H = [4  NaN  NaN]; %% Translocation rate [gC /m^2 d]
LtR_H = [1  NaN  NaN]; %%% Leaf to Root ratio maximum
Mf_H= 1./[50  NaN  NaN]; %% fruit maturation turnover [1/d]
Wm_H= [0  NaN  NaN] ; % wood turnover coefficient [1/d]
eps_ac_H = [1  NaN  NaN]; %% Allocation to reserve parameter [0-1]
Klf_H = 1./[10  NaN  NaN]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74           NaN           NaN]; %% fraction above-ground sapwood and reserve
fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1           NaN           NaN]; %% Reference allocation to Fruit and reproduction
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
[Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
[Stoich_H(3)]=Veg_Stoichiometric_Parameter(Nl_H(3));
%%%% Phenology 
Bfac_lo_H= [0.96           NaN           NaN]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN  NaN  NaN] ;  % Not-used 
Tlo_H = [10  NaN  NaN]; %% Mean Temperature for Leaf onset
Tls_H = [NaN  NaN  NaN]; %%% Not-used 
PAR_th_H= [NaN          NaN           NaN]; %% Light Phenology Threshold 
dmg_H= [10  NaN  NaN]; %%%  Day of Max Growth
LAI_min_H = [0.1           NaN           NaN];
mjDay_H = [-1  NaN  NaN]; %% Maximum Julian day for leaf onset
LDay_min_H =[10.0           NaN           NaN]; %% Minimum Day duration for leaf onset
LDay_cr_H = [10  NaN  NaN]; %%%  Threshold for senescence day light [h]
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
Sl_L = [NaN         0.012         0.020]; %  [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_L = [NaN   32   45]; %[kgC/kgN ] Leaf Nitrogen Concentration
r_L =  [NaN         0.055          0.06];  %% respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_L=  [NaN          0.25          0.25];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
dd_max_L= 1./[NaN   10   12];%%%  [1/d]  death maximum for drought
dc_C_L =  1./[NaN      52          55]; %%[Factor of increasing mortality for cold]
Tcold_L = [NaN    2    0]; %% [°C] Cold Leaf Shed
drn_L=  1./[NaN   450  1290]; %% turnover root  [1/d]
dsn_L= 1./[NaN  365  365]; % normal transfer rate sapwood [1/d]
age_cr_L= [NaN  180  180]; %% [day] Critical Leaf Age
Trr_L = [NaN           5.0           5.0]; %% Translocation rate [gC /m^2 d]
LtR_L = [NaN           0.7           0.7]; %%% Leaf to Root ratio maximum
Mf_L= 1./[NaN   50   50]; %% fruit maturation turnover [1/d]
Wm_L= [NaN    0    0] ; % wood turnover coefficient [1/d]
eps_ac_L = [NaN           1.0         1.0]; %% Allocation to reserve parameter [0-1]
Klf_L = 1./[NaN   20   20]; %% Dead Leaves fall turnover [1/d]
fab_L = [NaN    0    0]; %% fraction above-ground sapwood and reserve
fbe_L = 1-fab_L; %% fraction below-ground sapwood and reserve
ff_r_L= [NaN           0.1           0.1]; %% Reference allocation to Fruit and reproduction
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
[Stoich_L(2)]=Veg_Stoichiometric_Parameter(Nl_L(2));
[Stoich_L(3)]=Veg_Stoichiometric_Parameter(Nl_L(3));
%%%% Phenology 
Bfac_lo_L= [NaN          0.99          0.99]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN  0.15 0.15] ; % Not-used 
Tlo_L = [NaN    8    6]; %% Mean Temperature for Leaf onset
Tls_L = [NaN  NaN  NaN]; %% Not-used 
PAR_th_L= [NaN  NaN    NaN]; %% Light Phenology Threshold 
dmg_L= [NaN   10   10]; %%%  Day of Max Growth
LAI_min_L = [NaN          0.02          0.01];
mjDay_L = [NaN   366  366]; %% Maximum Julian day for leaf onset
LDay_min_L =[NaN   10   10]; %% Minimum Day duration for leaf onset
LDay_cr_L = [NaN   10   10]; %%%  Threshold for senescence day light [h]
[ParEx_L(1)]=Exudation_Parameter(0);
[Mpar_L(1)]=Vegetation_Management_Parameter;
[ParEx_L(2)]=Exudation_Parameter(0);
[Mpar_L(2)]=Vegetation_Management_Parameter;
[ParEx_L(3)]=Exudation_Parameter(0);
[Mpar_L(3)]=Vegetation_Management_Parameter;
%%%%%%%%%%
% Mpar_L(2).jDay_cut=[1:365];
% Mpar_L(2).LAI_cut=[-0.1]; %% Grazing %%
% Mpar_L(3).jDay_cut=[1:365];
% Mpar_L(3).LAI_cut=[-0.1]; %% Grazing %%
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
LAI_H(1,:)=[0.0  0  0]; %
B_H(1,:,:)= [6 71 42 26 2 0 7 0 ;0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0]; %% 
Rrootl_H(1,:)= [1073  0  0] ;
PHE_S_H(1,:)=[3  0  0];
dflo_H(1,:)=[0  0  0];
AgeL_H(1,:)=[545  0  0];
e_rel_H(1,:)=[1  1  1];
hc_H(1,:) =[5.2  0  0]; %%
SAI_H(1,:) = [0.2  0  0]; %% 

Bfac_weekH(1,:)=[0 0 0]; 
Bfac_dayH(1,:)=[0 0 0]; 
%%%%%%%%%%%%%%%%%%%%%
LAI_L(1,:)=[0  0  0]; %
B_L(1,:,:)= [0 0 0 0 0 0 0 0 ; 0 0 135 137 6 0 1 0  ;0 0 56 44 1 0 1 0]; %% 
Rrootl_L(1,:)= [0  578 609] ;
PHE_S_L(1,:)=[0  1  1];
dflo_L(1,:)=[0  0  0];
AgeL_L(1,:)=[0  0  0];
e_rel_L(1,:)=[1  1  1];
hc_L(1,:) =[0  0.2  0.15]; %%
SAI_L(1,:) = [0  0.001  0.001]; %% 

Bfac_weekL(1,:)=[0 0 0]; 
Bfac_dayL(1,:)=[ 0 0 0];
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
O(1,:)=[   0.0411    0.0296    0.0161    0.0161    0.0161    0.0161    0.0161    0.0161    0.0161    0.0476    0.1184    0.1505    0.1606    0.1658,...  
    0.1695    0.1725    0.1748  0.1773]; 
%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%%%%%%%%%%%%%%%%%
cd(cur_dir)
%%%%%%%%%%%%%%%%%