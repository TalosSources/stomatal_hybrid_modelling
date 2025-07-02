%%%%%%%%%%%%%%%%%%%%%%% Wallaby %%%%%%%%%%%%%
%%%%%%%%%%% Pr = 1200 (1080) Ta = 12-22  Ws 2.4-2.6
% Tall wet sclerophyll forest, Eucalyptus regnans dominant,
%canopy height 80–90 m, prominent rainforest
%(dense) understorey wet sclerophyll with Pomaderris aspera and Olearia argophylla, which
%are typically 10–18 m tall, and tree ferns (Cyathea australis and Dicksonia antartica)
%moss-covered logs.
%%% ET 814-911 mm/yr GPP 1980 (2615 pre fire)
% E. regnans only 36.8 stems per hectare, there are 37±2 stems ha?1 of Mountain Ash, with a mean
%diameter at breast height (DBH) of 199±8 cm and canopy height of 80±5 m, though some individual trees are taller than 90 m (van Pelt et al., 2004).
%Below the Mountain Ash open canopy, there is a dense temperate rainforest understorey (298±14 stems ha?1)

%%% Litterfall 60-100 gDM/m2/month  (400 gDM/m2 yr) 4 ton/ha/yr   %%
%%% Mean litter and CWD were estimated to be 5 and 26 t C ha?1 respectively
%%% Fine roots 650-720 gC/m2  Ra-roots 1070 gC/m2
%%% LAI = 1.69 - 2.5  eucalyptos  -- 3.9 the LAI for the site could be as
%%% high 10.0 m2m?2  1.7 + 1.7
%Aboveground biomass (AGB) is between 700 and 920 t C ha?1


%% Eucalyptus regnans;; Temperate Rainforest species
% The soils of Kinglake and the Hume plateau are rich krasnozemic
%loams, which are friable, red brown gradational soils with high
%organic matter contents (15--20%) in the upper 20--30 cm (Ashton
%2000b). The soil profile at 0--30 cm is a dark brown coarsely
%friable loam and grades with depth to a red-brown or yellowbrown
%subsoil clay loam. From a depth of 183--244 cm the soil
%becomes a sandy loam and from 244--274 cm the soil becomes
%loamy gravel.

%%%%%% SOIL AND HYDROLOGICAL PARAMETER
%%%%%%%%%%%%%%%%%
cur_dir=cd;
cd(Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VEG. SPECIES
%%% Rainfall Disaggregation
a_dis = NaN ;
pow_dis = NaN;
%%%%%%%%%%%%%%%%%%%%
fpr=1;
SvF=1; %% Sky View Factor
SN=0; %% Stream Identifier
Slo_top=0;  %% [fraction dy/dx]
Slo_pot=zeros(1,ms); %% [fraction dy/dx]
Asur = 1./cos(atan(Slo_top)); %% Real Area/Projected Area [m^2/m^2]
Ared = 1;
aR =1; %%% anisotropy ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsize=1; %%[m^2];
aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
Kbot =  NaN; %% [mm/h] Conductivity at the bedrock layer
Krock = NaN; %% [mm/h] Conductivity of Fractured Rock
zatm = 100; %% Reference Height
%%% zatm = 5; %% After Fire ;;
%%%%%%%%%%%%%%%%%%
%%%% LAND COVER PARTITION
%%% Eucalyptos Regnans // Understory
Cwat = 0; Curb = 0.0 ; Crock = 0.0;
Cbare = 0; Ccrown = [0.65 0.35];
%%%%%%%%%% SOIL INPUT %%% Sandy soil
Pcla= 0.32;
Psan= 0.30;
Porg= 0.025;
Color_Class = 0;
%%%%%%%%%%%%%%%%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%%%%%%%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s =lan_s*ones(1,ms);
cv_s = cv_s*ones(1,ms);
%%%%
Ks=Ks+10;
%SPAR=2;
SPAR=2; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
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
%%%%
%nVG(7:end)=1.28;
%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
clear Oss_Lp Owp_Lp
Oice = 0;
%%%
Zs= [0 10 50 100 150 200 300 400 500 600 700 800 900 1000 1100 1300 1500 2000 2500 3000 3500 4000 4500 ]; %%% [ms+1]
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
%%%%%%%%%%
Ks_Zs=Ks*exp(-cumsum(Dz)/1000); 
%%%%%%%%%%%%%%%%% OTHER PARAMETER
In_max_urb=5; In_max_rock=0.1; %% [mm]
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExEM = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS VEGETATION
%%% Root Depth
%%% cc -- number of crown area
CASE_ROOT= 1 ;  %%% Type of Root Profile
ZR95_H = [3500 2200]; %% [mm]
ZR95_L = [0 0]; %% [mm]
ZR50_H = [NaN NaN];
ZR50_L =[NaN NaN];
ZRmax_H = [NaN NaN];
ZRmax_L = [NaN NaN];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_L_In= [0.2 0.2]; %%[mm/LAI]
Sp_LAI_H_In= [0.2 0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [5 2]; %%[cm]
d_leaf_L= [0 0];  %% [cm]
%%%%%%%% Biochemical parameter
KnitH= [0.25 0.30]; %%% Canopy Nitrogen Decay
KnitL= [0.0 0.0];
mSl_H = [0.0 0.0];%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = [0.0 0.0]; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Photosynthesis Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_H=[0.081 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[1000 1000]; %%[Pa]
a1_H=[8 9];% [8 8];
go_H=[0.01 0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3 3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[45 45] ;% [72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[2.0 2.0]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FI_L=[0.040 0.040];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[1800 1800]; %%[Pa]
a1_L=[6 6];
go_L=[0.01 0.01];% % [mol / s m^2] minimum Stomatal Conductance
CT_L=[4 4];  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[0.649 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[Inf Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L=[2.1 2.1]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pss_H = [1000]; Pwp_H = [3800]; %%% [kPa]
%Pss_L = [850]; Pwp_L = [3000]; %%% [kPa]
%%%%%%%% Hydraulic Parameters
%%% Stomata -1 -3.8
Psi_sto_00_H =  [-0.7 -0.8]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H =  [-1.5 -2.0];%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H = [ -1.4 -1.0]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H = [ -2.5 -2.5];%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H = [10 10] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200 1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = [10 10.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [ -9.0 -9.0]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [150 150]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata
Psi_sto_00_L =  [-0.5 -0.5]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_L =  [-1.6 -1.6] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L =  [-1.0 -1.0]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_L =  [-2.0 -2.0] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [ 5 5] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [1200 1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L = [6.0 6.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [-10 -10.0]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [150 150]; %%% [kg / m^3 sapwood MPa]
%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);

%%%% Growth Parameters
PsiG50_H= [-0.7 -0.8];  %%[MPa]
PsiG99_H= [-1.5 -2.0];  %%[MPa]
gcoef_H = [3.5 3.5]; % [gC/m2 day]
%%------
PsiG50_L= [-1.6;-1.6];  %%[MPa]
PsiG99_L= [-2.0; -2.0] ; %%[MPa]
gcoef_L = [3.5 3.5]; % [gC/m2 day]

%%%%%%%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(5);
[PFT_opt_H(2)]=Veg_Optical_Parameter(5);
[PFT_opt_L(1)]=Veg_Optical_Parameter(0);
[PFT_opt_L(2)]=Veg_Optical_Parameter(0);
OM_H=[1 1];
OM_L=[1 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.012 0.018]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [60 42]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
[Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
r_H = [0.042 0.055];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_H= [0.25 0.25]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_H= [2.0 10.0];
aSE_H= [0 0]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_H= [1/365 1/150]; %%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  [78/365 78/365]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H = [-15 -15]; %% [°C] Cold Leaf Shed
drn_H=  [1/600 1/700]; %% turnover root  [1/d]
dsn_H= [1/450 1/550]; % normal transfer rate sapwood [1/d]
age_cr_H= [365 365]; %% [day] Critical Leaf Age
Bfac_lo_H= [0.98 0.99]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [8.8 8.5]; %% Mean Temperature for Leaf onset
Tls_H = [NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_H= [NaN NaN];
dmg_H= [20 25]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.001 0.001];
Trr_H = [1.0 0.5]; %% Translocation rate [gC /m^2 d]
mjDay_H = [-240 -240]; %% Maximum Julian day for leaf onset
LDay_min_H =[11.4 11.4]; %% Minimum Day duration for leaf onset
LtR_H = [0.8 0.6]; %%% Leaf to Root ratio maximum
Mf_H= [1/60 1/60]; %% fruit maturation turnover [1/d]
Wm_H= [0 0] ; % wood turnover coefficient [1/d]
eps_ac_H = [0.6 0.5]; %% Allocation to reserve parameter [0-1]
LDay_cr_H = [10.6 10.6]; %%%  Threshold for senescence day light [h]
Klf_H =[1/40 1/40]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74 0.74]; %% fraction above-ground sapwood and reserve
fbe_H = [0.26 0.26]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1 0.1]; %% Reference allocation to Fruit and reproduction
[ParEx_H(1)]=Exudation_Parameter(0);
[ParEx_H(2)]=Exudation_Parameter(0);
[Mpar_H(1)]=Vegetation_Management_Parameter;
[Mpar_H(2)]=Vegetation_Management_Parameter;
% %%%%
% Mpar_H(1).Date_fire =  datenum(2009,2,7,1,0,0) ; %%% Bushfire 7 February 2009
% Mpar_H(1).fire_eff = 1;
% Mpar_H(1).fract_left=0;
% Mpar_H(1).funb_nit= 0.25; %%
% Mpar_H(1).fract_resprout =0.0;% 1.0;
% %%%%
% Mpar_H(2).Date_fire =  datenum(2009,2,7,1,0,0) ; %%% Bushfire 7 February 2009
% Mpar_H(2).fire_eff = 1;
% Mpar_H(2).fract_left=0;
% Mpar_H(2).funb_nit= 0.25; %%
% Mpar_H(2).fract_resprout = 0.0;% 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L = [0.020 0.020]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L= [28 28]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
[Stoich_L(2)]=Veg_Stoichiometric_Parameter(Nl_L(2));
%PLNR_L = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_L = [0.050 0.050];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_L= [0.25 0.25]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_L= [2.0 10.0];
aSE_L= [2 2]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_L= [1/12 1/12]; %[1/50 1/50];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L = [7/365 7/365]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_L = [2.0 2.0]; %% [°C] Cold Leaf SLed
drn_L= [1/650 1/650 ]; %% [1/900 1/1200]; turnover root  [1/d]
dsn_L= [1/365  1/365 ]; % normal transfer rate sapwood [1/d]
age_cr_L= [150 150]; %% [day] Critical Leaf Age
Bfac_lo_L= [0.95 0.95]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_L = [8.0 8.0]; %% Mean Temperature for Leaf onset
Tls_L = [NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_L= [-Inf -Inf];
dmg_L= [20 20]; %%% 20 -- Tree 30 Grasses Day of Max Growth
LAI_min_L = [0.05 0.05 ];
Trr_L = [3.5 3.5]; %% 1.3 Translocation rate [gC /m^2 d]
mjDay_L = [-1 -1]; %% Minimum Julian day for leaf onset
LDay_min_L =[10.2 10.0]; %% Day duration for leaf onset
LtR_L =  [0.6 0.6]; %%% Leaf to Root ratio maximum
Mf_L= [1/50 1/50]; %% fruit maturation turnover [1/d]
Wm_L= [0 0] ; % wood turnover coefficient [1/d]
eps_ac_L = [1.0 1.0]; %% [0.5] Allocation to reserve parameter [0-1]
LDay_cr_L = [10.6 10.6] ; %[11.0 11.0]; %%%  Threshold for senescence day light [h] [12.4 12.4]
Klf_L = [1/20 1/20]; % [1/100 1/100];% ; %% Dead Leaves fall turnover [1/d]
fab_L = [0.0 0.0 ]; %% fraction above-ground sapwood and reserve
fbe_L = [1.0 1.0 ]; %% fraction below-ground sapwood and reserve
ff_r_L= [0.1 0.1];
[ParEx_L(1)]=Exudation_Parameter(0);
[ParEx_L(2)]=Exudation_Parameter(0);
[Mpar_L(1)]=Vegetation_Management_Parameter;
[Mpar_L(2)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%% PRODUCTIVITY
%%% Maximum Rubisco Capacity
%[Vmax_H]=Maximum_Rubisco_Capacity(Sl_H,PLNR_H,Nl_H);
%[Vmax_L]=Maximum_Rubisco_Capacity(Sl_L,PLNR_L,Nl_L);
Vmax_H = [50 46]; %55
Vmax_L = [0 0]; %35
%[Amax_H]= MAX_PHOTOSYNTESIS(Vmax_H,Ca,CT_H,Tup_H,Tlow_H,FI_H,Oa); %% 17
%[Amax_L]= MAX_PHOTOSYNTESIS(Vmax_L,Ca,CT_L,Tup_L,Tlow_L,FI_L,Oa); %% 13
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')

%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%%%%%%%%%%
LAI_H(1,:)=[4.49 4.71]; %
B_H(1,:,:)= [374 706 467 473 24 0 39 0 ; 274 734 458 491 21 0 29 0]; %%
%B_H(1,:,:)= [380 660 475 445 22 31400 40 0 ; 290 690 485 463 21 7400 30 0]; %%
Rrootl_H(1,:)= [4907 5246] ;
PHE_S_H(1,:)=[3 3];
dflo_H(1,:)=[29 33];
AgeL_H(1,:)=[344 343];
e_rel_H(1,:)=[1 1];
hc_H(1,:) =[80 18]; %%
SAI_H(1,:) = [0.2 0.2]; %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LAI_L(1,:)=[0.0 0.0];
B_L(1,:,:)= [0 0 0 0 0 0 0 0 ; 0 0 0 0 0 0 0 0 ];
Rrootl_L(1,:)= [0 0] ;
PHE_S_L(1,:)=[3 3];
dflo_L(1,:)=[0 27 ];
AgeL_L(1,:)=[0 16];
e_rel_L(1,:)=[1 1];
hc_L(1,:) =[0.0 0.0];
SAI_L(1,:) = [0.0 0.0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [1000 1000];
Preserve_H(1,:)= [1000 1000];
Kreserve_H(1,:)= [1000 1000];
FNC_H(1,:)=[1 1];
NupI_H(1,:,:)= [0 0 0 ; 0 0 0];
Nreserve_L(1,:)= [1000 1000];
Preserve_L(1,:)= [100 1000];
Kreserve_L(1,:)= [100 1000];
FNC_L(1,:)=[1 1];
NupI_L(1,:,:)= [0 0 0 ; 0 0 0];
RexmyI(1,:)= [0 0 0];
%%%
%%%%%
if OPT_SoilBiogeochemistry == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Wm_H= [1/33580 1/21900];
    B_H(1,:,:)= [374 706 467 473 24 42432 39 0 ; 274 734 458 491 21 12000 29 0]; %%
    Nreserve_H(1,:)= [8.0 10.1];
    Preserve_H(1,:)= [0.22 0.22];
    Kreserve_H(1,:)= [1.0 1.0];
    FNC_H(1,:)=[1 1];
    NupI_H(1,:,:)= [  0.0189 0.00119  0.0165 ;
        0.0334 0.00222   0.0218];
    RexmyI(1,:)= [   0.047149422973637   0.141208060799000      0];
    Nreserve_L(1,:)= [0 0]; Preserve_L(1,:)= [0 0]; Kreserve_L(1,:)= [0 0]; FNC_L(1,:)=[1 1];
    NavlI(1,:)=[   0.419441070947630   2.238857286638924   0.163542877995922];
    
    %%%%%%
   P(1,:)=  10000*[   0.007275634661199   0.038856949512524   0.006381202247148   0.259132253643716   0.086377417881238   0.005314894022662,...
   0.051310795551702   0.010419073104819   0.436324912510113   0.150353859404462   1.075761899519711   0.001489499191890,...
   0.001202316833399   0.000010357797413   0.000006648758034   0.000004315748922   0.000011081263389   0.005254156411233,...
   0.017992615705609   0.005560393781658  -0.000000000000557   0.000406545445558   0.000238407583276   0.000915180313728,...
   0.000268772845359   0.101634095835779   0.001010667435953   0.002770082093153   0.000308442789262  -0.000000000000274,...
   0.000003457767642   0.000000357643224   0.000002031960699   0.000040654544556   0.000041992356369   0.000065370022409,...
   0.000031810694506   0.021965251828231   0.000327925985640   0.000450526628085   0.000046266418389  -0.000000000000274,...
   0.000000808984956   0.015000000000000   0.000051019086837   0.001500269183892   0.000000028087167   0.000388716652698,...
   0.000457590156864   0.000075718933939   0.001626089892707   0.000043099216775   0.000042307531059   0.000728830276879,...
   0.050325675876242];


    
    Upl=0.01; %% Soil production [mm/yr]
    HIST=0;
    %[B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl,HIST,FertN,FertP,FertK,DepN,DepP,DepK);
    [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl,HIST);
    B_IO.SC_par=[1 1 1 1];
    PHs=5.9;  %
end
%%%%
TBio_Lt(1,:)=[0 0];  %%[ton DM / ha ]
TBio_Ht(1,:)=[1800 10]; %[ton DM / ha ]
AgePl_H(1,:)= [300*365.25 20*365.25]; %% days
%%%%%%%%%%%%%%%%%
Vx_H=[10 10];  %% [mm/ m2 PFT];
Vl_H=[10 10];  %% [mm/ m2 PFT];
Vx_L=[10 10];   %% [mm/ m2 PFT];
Vl_L=[10 10];   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%%%%%%% Initial Conditions
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+2;
Tdamp(1)=14;
Tdp(1,:)= Ta(1)*ones(1,ms);
TdpI_H(1,:)=14.0;
TdpI_L(1,:)=14.0;
%%% Snow_alb = soil_alb initial
snow_alb.dir_vis = 0.2;
snow_alb.dif_vis = 0.2;
snow_alb.dir_nir = 0.2;
snow_alb.dif_nir = 0.2;
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
Ws_under(1,:)=1;
%%%%%%%%%%%%%% Volume [mm]
%%%%%%%%%%%%%%%%%%%%%
O(1,:)= [ 0.3043    0.3045    0.3049    0.3054    0.3060    0.3070    0.3082    0.3095    0.3107    0.3118    0.3129    0.3140,...
    0.3150    0.3160    0.3175    0.3193    0.3224    0.3263    0.3295   0.3321    0.3340    0.3349];

%dz(ms)=ZWT(1)-Zs(ms);
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%%%%%%%%%%%%%%%%%
cd(cur_dir)
%%%%%%%%%%%%%%%%%