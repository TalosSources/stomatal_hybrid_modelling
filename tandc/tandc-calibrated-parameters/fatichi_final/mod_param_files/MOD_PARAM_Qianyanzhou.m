% %%%% Qianyanzhou planted coniferous forest
%%% planted in 1983.The dominant species are Pinus elliottii,Pinus massoniana,Cunninghamia lanceolata and Schima superba
%%% Pr is 1488.8mm,and evaporation is 1110.3 mm
%  A 42 m-height flux observation tower / /Canopy height 13m 
%%% typical red soil. Sand 27% Clay 19% 
%the leaf area index (LAI) of the plantation was 4.5  (3.6-4.9 seasonality)
% dBiomass = 450 gC m-2 yr-1  GPP = 1810 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETER FILE TEMPLATE FOR T&C %%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
cur_dir=cd;
cd(Directory)
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
zatm =42; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%
%%%% LAND COVER and VEGETATION COMPOSITION 
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare = 0.0; Ccrown = [1];
%%% Rainfall disaggrgation information 
a_dis = NaN ;
pow_dis = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOIL PARAMETERS 
Pcla= 0.19;
Psan= 0.27;
Porg= 0.02;
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
Zs= [0    10    20    50   100   150   200   300   400   600  700  800  900 1000  1250  1500 ]; %% ms+1
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
TminS=-0.8;%% Threshold temperature snow
TmaxS=2.8;%% Threshold temperature snow
ros_max1=580; %%% [kg/m^3]
ros_max2=300; %%% [kg/m^3]
Th_Pr_sno = 8; %%% [mm/h] Threshold Intensity of snow to consider a New SnowFall
%%%% ICE Parameter
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
Aice = 0.28; %% [-] Ice albedo
WatFreez_Th = -8; %% [°C] Threshold for freezing water
dz_ice = 0.45; %% [mm/h] Freezing Layer depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%  ROOT PARAMETER 
ExEM = 1; %% Fraction of ectomychorrizae per area 
CASE_ROOT= 1;  %%% Type of Root Profile
ZR95_H = [1200]; %% [mm]
ZR95_L = [0]; %% [mm]
ZR50_H = [NaN];
ZR50_L = [NaN];
ZRmax_H = [NaN];
ZRmax_L = [NaN];
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
Sp_LAI_H_In= [0.1]; %%[mm/LAI]
Sp_LAI_L_In= [0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [0.25]; %%[cm]
d_leaf_L= [2];  %% [cm]
%% Veg Biochemical parameter
KnitH=[0.40]; %%% Canopy Nitrogen Decay
KnitL=[NaN];
%%%%%
mSl_H = [0];%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = [NaN]; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%
%%------
%%------
FI_H=[0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[800]; %%[Pa] 
a1_H=[5];  %%% [-] WUE parameter 
go_H=[0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[2.0]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------
FI_L=[NaN];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[NaN]; %%[Pa] 
a1_L=[NaN];  %%% [-] WUE parameter 
go_L=[NaN];% [mol / s m^2] minimum Stomatal Conductance
CT_L=[NaN]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[NaN];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[NaN]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[NaN]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L=[NaN]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%
Vmax_H = [45]; % [umol CO2 /m2 s] - Maximum Rubisco Capacity 
Vmax_L = [0]; % [umol CO2 /m2 s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hydraulic Parameters
Psi_sto_00_H = [-0.8]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H = [-1.8] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H = [-1.5]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H = [-3.2] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H = [20] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_H = [6] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-8.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [80]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata
Psi_sto_00_L= [NaN]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_L = [NaN] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L = [NaN]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_L = [NaN] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [NaN] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [NaN];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_L = [NaN] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [NaN];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [NaN]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [NaN]; %%% [kg / m^3 sapwood MPa]
%% Growth Parameters
PsiG50_H= [-0.8];  %%[MPa]
PsiG99_H= [-1.8];  %%[MPa]
gcoef_H = [3.5]; % [gC/m2 day]
%%------
PsiG50_L= [NaN];
PsiG99_L= [NaN];
gcoef_L = [NaN]; % [gC/m2 day]
%%%
%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(1);
[PFT_opt_L(1)]=Veg_Optical_Parameter(0);
OM_H=[1];
OM_L=[NaN];
%%%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  VEGETATION PART %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% High Vegetation 
aSE_H=  [0]; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_H = [0.018]; %  [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_H = [42]; %[kgC/kgN ] Leaf Nitrogen Concentration
r_H =  [0.055];  %% respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_H=  [0.25];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
dd_max_H= 1./[200];%%%  [1/d]  death maximum for drought
dc_C_H =  1./[122]; %%[Factor of increasing mortality for cold]
Tcold_H = [2]; %% [°C] Cold Leaf Shed
drn_H=  1./[500]; %% turnover root  [1/d]
dsn_H= 1./[450]; % normal transfer rate sapwood [1/d]
age_cr_H= [550]; %% [day] Critical Leaf Age
Trr_H = [0.7]; %% Translocation rate [gC /m^2 d]
LtR_H = [0.9]; %%% Leaf to Root ratio maximum
Mf_H= 1./[50]; %% fruit maturation turnover [1/d]
Wm_H= [0] ; % wood turnover coefficient [1/d]
eps_ac_H = [0.2]; %% Allocation to reserve parameter [0-1]
Klf_H = 1./[25]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74]; %% fraction above-ground sapwood and reserve
fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1]; %% Reference allocation to Fruit and reproduction
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
%%%% Phenology 
Bfac_lo_H= [0.99]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN] ;  % Not-used 
Tlo_H = [15]; %% Mean Temperature for Leaf onset
Tls_H = [NaN]; %%% Not-used 
PAR_th_H= [NaN]; %% Light Phenology Threshold 
dmg_H= [45]; %%%  Day of Max Growth
LAI_min_H = [0.001];
mjDay_H = [170]; %% Maximum Julian day for leaf onset
LDay_min_H =[12.2]; %% Minimum Day duration for leaf onset
LDay_cr_H = [12.2]; %%%  Threshold for senescence day light [h]
[ParEx_H(1)]=Exudation_Parameter(0);
[Mpar_H(1)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Low Vegetation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
aSE_L=  [NaN]; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_L = [NaN]; %  [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_L = [NaN]; %[kgC/kgN ] Leaf Nitrogen Concentration
r_L =  [NaN];  %% respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_L=  [NaN];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
dd_max_L= 1./[NaN];%%%  [1/d]  death maximum for drought
dc_C_L =  1./[NaN]; %%[Factor of increasing mortality for cold]
Tcold_L = [NaN]; %% [°C] Cold Leaf Shed
drn_L=  1./[NaN]; %% turnover root  [1/d]
dsn_L= 1./[NaN]; % normal transfer rate sapwood [1/d]
age_cr_L= [NaN]; %% [day] Critical Leaf Age
Trr_L = [NaN]; %% Translocation rate [gC /m^2 d]
LtR_L = [NaN]; %%% Leaf to Root ratio maximum
Mf_L= 1./[NaN]; %% fruit maturation turnover [1/d]
Wm_L= [NaN] ; % wood turnover coefficient [1/d]
eps_ac_L = [NaN]; %% Allocation to reserve parameter [0-1]
Klf_L = 1./[NaN]; %% Dead Leaves fall turnover [1/d]
fab_L = [NaN]; %% fraction above-ground sapwood and reserve
fbe_L = 1-fab_L; %% fraction below-ground sapwood and reserve
ff_r_L= [NaN]; %% Reference allocation to Fruit and reproduction
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
%%%% Phenology 
Bfac_lo_L= [NaN]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN] ; % Not-used 
Tlo_L = [NaN]; %% Mean Temperature for Leaf onset
Tls_L = [NaN]; %% Not-used 
PAR_th_L= [NaN]; %% Light Phenology Threshold 
dmg_L= [NaN]; %%%  Day of Max Growth
LAI_min_L = [NaN];
mjDay_L = [NaN]; %% Maximum Julian day for leaf onset
LDay_min_L =[NaN]; %% Minimum Day duration for leaf onset
LDay_cr_L = [NaN]; %%%  Threshold for senescence day light [h]
[ParEx_L(1)]=Exudation_Parameter(0);
[Mpar_L(1)]=Vegetation_Management_Parameter;
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
LAI_H(1,:)=[4.38]; %
B_H(1,:,:)= [243 489 270 336 1 0 10 0]; %% 
Rrootl_H(1,:)= [0] ;
PHE_S_H(1,:)=[1];
dflo_H(1,:)=[0];
AgeL_H(1,:)=[537];
e_rel_H(1,:)=[1];
hc_H(1,:) =[13]; %%
SAI_H(1,:) = [0.1]; %% 
%%%%%%%%%%%%%%%%%%%%%
LAI_L(1,:)=[0]; %
B_L(1,:,:)= [0 0 0 0 0 0 0 0]; %% 
Rrootl_L(1,:)= [0] ;
PHE_S_L(1,:)=[0];
dflo_L(1,:)=[0];
AgeL_L(1,:)=[0];
e_rel_L(1,:)=[1];
hc_L(1,:) =[0]; %%
SAI_L(1,:) = [0]; %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [1000];
Preserve_H(1,:)= [1000];
Kreserve_H(1,:)= [1000];
FNC_H(1,:)=[1];
NupI_H(1,:)= [0 0 0];
Nreserve_L(1,:)= [1000];
Preserve_L(1,:)= [1000];
Kreserve_L(1,:)= [1000];
FNC_L(1,:)=[1];
NupI_L(1,:)= [0 0 0];
RexmyI(1,:)= [0 0 0];
%%%
%%%%
if OPT_SoilBiogeochemistry == 1
    %%%%%
end
%%%
TBio_L=[1];  %%[ton DM / ha ]
TBio_H=[1]; %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=[1];  %% [mm/ m2 PFT];
Vl_H=[1];  %% [mm/ m2 PFT];
Vx_L=[1];   %% [mm/ m2 PFT];
Vl_L=[1];   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%%%% Initialization OTHER
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+1;
Tdamp(1)=Ta(1);
Tdp(1,:)= 11*ones(1,ms);
TdpI_H(1,:)=11;
TdpI_L(1,:)=11;
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
O(1,:)= Ofc; 
%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%%%%%%%%%%%%%%%%%
cd(cur_dir)
%%%%%%%%%%%%%%%%%