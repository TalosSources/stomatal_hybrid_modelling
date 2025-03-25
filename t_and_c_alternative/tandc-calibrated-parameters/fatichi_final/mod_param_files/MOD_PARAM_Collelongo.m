%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETER FILE TEMPLATE FOR T&C %%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Site:	Collelongo				
% beech (Fagus sylvatica L.) forest
% It is a pure beech forest with 830 trees ha?1, with an average
% diameter of 22 cm and an average height of 21.5 m. The average tree age is approximately
% 120 years
% 
% Tree species:	beech (Fagus sylvatica L.)	
% Canopy height:	22 m (average); 25 m (highest trees)	
% Canopy foliage:	Generally 80% of LAI between 13 and 22 m	
% Eddy covariance sampling height:	32 m
%% GPP=1258 - 1345   AG - Biomass 13700 gC/m2   LAI = 5.0 
%Sand Silt Clay  (%) 6 50 44  
%Root:shoot 0.28  Literfall = 245 gC/m2   AG-Biomass growth 363 gC/m2 
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
zatm =32; %% Reference Height
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
Pcla= 0.44;
Psan= 0.06;
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
%%%%%%%%%%%%%
Zs= [0    10    20    50   100   150   200   300   400   500   600   800  1000  1250  1500  1750  2000]; %% ms+1
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
TminS=-0.7;%% Threshold temperature snow
TmaxS=2.8;%% Threshold temperature snow
ros_max1=580; %%% [kg/m^3]
ros_max2=300; %%% [kg/m^3]
Th_Pr_sno = 8; %%% [mm/h] Threshold Intensity of snow to consider a New SnowFall
%%%% ICE Parameter
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
Aice = 0.28; %% [-] Ice albedo
WatFreez_Th = -8; %% [°C] Threshold for freezing water
dz_ice = 0.54; %% [mm/h] Freezing Layer depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%  ROOT PARAMETER 
ExEM = 0; %% Fraction of ectomychorrizae per area 
CASE_ROOT= 1;  %%% Type of Root Profile
ZR95_H = [1500]; %% [mm]
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
Sp_LAI_H_In= [0.15]; %%[mm/LAI]
Sp_LAI_L_In= [0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [4]; %%[cm]
d_leaf_L= [2];  %% [cm]
%% Veg Biochemical parameter
KnitH=[0.25]; %%% Canopy Nitrogen Decay
KnitL=[NaN];
%%%%%
mSl_H = [0.0];%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = [NaN]; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%
%%------
%%------
FI_H=[0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[1000]; %%[Pa] 
a1_H=[6];  %%% [-] WUE parameter 
go_H=[0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[76]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[2.4]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
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
Vmax_H = [62]; % [umol CO2 /m2 s] - Maximum Rubisco Capacity 
Vmax_L = [0]; % [umol CO2 /m2 s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hydraulic Parameters
Psi_sto_00_H = [-0.8]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H = [-2.5] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H = [-1.2]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H = [-3.5] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H = [10] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_H = [15] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-5.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [150]; %%% [kg / m^3 sapwood MPa]
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
PsiG99_H= [-2.5];  %%[MPa]
gcoef_H = [4.5]; % [gC/m2 day]
%%------
PsiG50_L= [NaN];
PsiG99_L= [NaN];
gcoef_L = [NaN]; % [gC/m2 day]
%%%
%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(7);
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
aSE_H=  [1]; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_H = [0.017]; %  [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_H = [28]; %[kgC/kgN ] Leaf Nitrogen Concentration
r_H =  [0.032];  %% respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_H=  [0.25];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
dd_max_H= 1./[100];%%%  [1/d]  death maximum for drought
dc_C_H =  1./[11]; %%[Factor of increasing mortality for cold]
Tcold_H = [0.0]; %% [°C] Cold Leaf Shed
drn_H=  1./[1200]; %% turnover root  [1/d]
dsn_H= 1./[800]; % normal transfer rate sapwood [1/d]
age_cr_H= [150]; %% [day] Critical Leaf Age
Trr_H = [5]; %% Translocation rate [gC /m^2 d]
LtR_H = [0.5]; %%% Leaf to Root ratio maximum
Mf_H= 1./[80]; %% fruit maturation turnover [1/d]
Wm_H= [0] ; % wood turnover coefficient [1/d]
eps_ac_H = [1]; %% Allocation to reserve parameter [0-1]
Klf_H = 1./[28]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74]; %% fraction above-ground sapwood and reserve
fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1]; %% Reference allocation to Fruit and reproduction
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
%%%% Phenology 
Bfac_lo_H= [0.95]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN] ;  % Not-used 
Tlo_H = [2]; %% Mean Temperature for Leaf onset
Tls_H = [NaN]; %%% Not-used 
PAR_th_H= [NaN]; %% Light Phenology Threshold 
dmg_H= [30]; %%%  Day of Max Growth
LAI_min_H = [0.01];
mjDay_H = [250]; %% Maximum Julian day for leaf onset
LDay_min_H =[12.8]; %% Minimum Day duration for leaf onset
LDay_cr_H = [10.8]; %%%  Threshold for senescence day light [h]
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
LAI_H(1,:)=[0];
B_H(1,:,:)= [ 0 526 444 366 14 0 72 0];
Rrootl_H(1,:)= [4500] ;
PHE_S_H(1,:)=[1];
dflo_H(1,:)=[0];
AgeL_H(1,:)=[0];
e_rel_H(1,:)=[1];
hc_H(1,:) =[22.0];
SAI_H(1,:) = [0.2];
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
Tdp(1,:)= Ta(1)*ones(1,ms);
Sdp(1,:)= (0.001*dz).*cv_s.*(Tdp(1,:)); %% [J °C/m^2 K]
TdpI_H(1,:)=2;
TdpI_L(1,:)=2;
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
