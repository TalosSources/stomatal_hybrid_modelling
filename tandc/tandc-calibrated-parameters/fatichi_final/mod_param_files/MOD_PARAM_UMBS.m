%%%%%%%%%%%% VEGETATION COMPOSITION
%Dominant canopy species include Populus grandidentata Michx.
%(bigtooth aspen) and Populus tremuloides Michx. (trembling
%aspen), which together comprise over 40% of the basal area,
%Quercus rubra L. (northern red oak), Betula papyrifera Marsh.
%(paper birch), Fagus grandifolia Ehrh. (American beech), Acer
%saccharum Marsh. (sugar maple), Acer rubrum L. (red maple), and
%Pinus strobus L. (white pine). The understory primarily consists
%of Pteridium aqulinium (bracken fern) and saplings of the canopy
%species.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cur_dir=cd;
cd(Directory)
%%% SOIL AND HYDROLOGICAL PARAMETER
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
%Kh=Ks*aR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsize=1; %%[m^2];
aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
Kbot = NaN; %% [mm/h] Conductivity at the bedrock layer
Krock = NaN; %% [mm/h] Conductivity of Fractured Rock
zatm = 46; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%%%%%% VEG. SPECIES  --- Low Grasses -- High Decidous
%%%% LAND COVER PARTITION
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare = 0.0; Ccrown = [1];
%%%%%%%%%% SOIL INPUT
Pcla= 0.006;
Psan= 0.926;
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
SPAR=1; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%% Measured %%%%
%Ks = 250;  Osat = 0.37;  nVG=2; alpVG= -0.004; lVG=0.5; %% Ohy=0.03;
Ks = 360;  Osat = 0.37;  nVG=1.68; alpVG= -0.0052; lVG=0.5; %% Ohy=0.03;
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
lVG= lVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%
nVG(7:end)=1.28;
%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
[s_SVG,bVG]=Soil_parameters_VG(Phy,Osat,Ohy,nVG,alpVG,0);
clear Oss_Lp Owp_Lp
Oice = 0;
%%%%%%%%%%%%%
Zs= [0 10 20 50 100 200 300 400 600 800 1000 1200 1500 1750 2000 2250 2500 3000 ]; %% ms+1
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
TminS=-0.7;%% Threshold temperature snow
TmaxS= 2.8;%% Threshold temperature snow
ros_max1=580; %600; %%% [kg/m^3]
ros_max2=300; %450; %%% [kg/m^3]
Th_Pr_sno = 8.0; %%% [mm/h] Threshold Intensity of snow to consider a New SnowFall
%%%% ICE Parameter
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
Aice = 0.28; %% [-] Ice albedo
WatFreez_Th = -8; %% [°C] Threshold for freezing water
dz_ice = 0.54; %% [mm/h] Freezing Layer depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
ExEM = 0.85;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS VEGETATION
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile
ZR95_H = [2000]; %% [mm]
ZR95_L = [0]; %% [mm]
ZR50_H = NaN;
ZR50_L = NaN;
ZRmax_H = NaN;
ZRmax_L = NaN;
%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
%5 Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_L_In= [0.0]; %%[mm/LAI]
Sp_LAI_H_In= [0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [4]; %%[cm]
d_leaf_L= [1];  %% [cm]
%%%%%%%% Biochemical parameter
KnitH=0.35; %%% Canopy Nitrogen Decay
KnitL=0.5;
%%%%%
mSl_H = 0.001;%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = 0.0; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%
%%------
%%------
FI_H=[0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[1000]; %%[Pa]
a1_H=[8];
go_H=[ 0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[76]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[2.4]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_L=[0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[ 2000]; %%[Pa]
a1_L=[9];
go_L=[0.02];% % [mol / s m^2] minimum Stomatal Conductance
CT_L=[3];  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[Inf];
rjv_L= [2.2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Hydraulic Parameters
%%% Stomata -0.8 -3.5
Psi_sto_00_H =  -0.8; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H =  -2.5 ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H =  -1.0; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H =  -3.5 ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H = 10 ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = 1200;  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = 15.0 ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = 80000;  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = -5.0; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= 150; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata
Psi_sto_00_L = -0.8;%  %% [MPa]  Water Potential at PLCs loss conductivity
Psi_sto_50_L = -3.0;%  %% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L =  -1.1; %% [MPa]  Water Potential at PLCs% loss conductivity
PsiL50_L =  -4.0 ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = 5 ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = 1200;  %%%  [500 - 3000]%  [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L = 0.0 ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = 80000;  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = -4.5; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= 150; %%% [kg / m^3 sapwood MPa]
%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);

%%%% Growth Parameters
PsiG50_H= -0.8;  %%[MPa]
PsiG99_H= -2.5;  %%[MPa]
gcoef_H = 4.5; % [gC/m2 day]
%%------
PsiG50_L= -1.45;
PsiG99_L= -4.0;
gcoef_L = 3.5; % [gC/m2 day]
%%%%%%%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(7);
[PFT_opt_L(1)]=Veg_Optical_Parameter(0);
OM_H=1;
OM_L=1;
%%%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.023]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [30]; %[kgC/kgN ] Leaf Nitrogen Concentration
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
%PLNR_H = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_H =  [0.030];  %% [0.066 -0.011] respiration rate at 10° [gC/gN d ]
gR_H=  [0.25];% [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_H=  [ NaN];
aSE_H=  [ 1]; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
dd_max_H=  [ 1/100];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  [ 36/365]; %%20 [Factor of increasing mortality for cold]
Tcold_H =  [4.5]; %%5 [°C] Cold Leaf Shed
drn_H=  [1/1900]; %% turnover root  [1/d]
dsn_H= [1/1000]; % normal transfer rate sapwood [1/d]
age_cr_H= [120]; %% [day] Critical Leaf Age
Bfac_lo_H= [0.95]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [1.5]; %% Mean Temperature for Leaf onset
Tls_H = [NaN]; %% Mean Temperature for Leaf Shed
PAR_th_H= [NaN];
dmg_H= [30]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.01];
Trr_H = [3.5]; %% Translocation rate [gC /m^2 d]
mjDay_H = [250]; %% Maximum Julian day for leaf onset
LDay_min_H =[13.90]; %% Minimum Day duration for leaf onset
LtR_H = [0.26]; %%% Leaf to Root ratio maximum
Mf_H= [1/80]; %% fruit maturation turnover [1/d]
Wm_H= [0] ; % wood turnover coefficient [1/d]
eps_ac_H = [1]; %% Allocation to reserve parameter [0-1]
LDay_cr_H = [11.5]; %%%  Threshold for senescence day light [h]
Klf_H =[1/28]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.80]; %% fraction above-ground sapwood and reserve
fbe_H = [0.20]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1]; %% Reference allocation to Fruit and reproduction
[ParEx_H(1)]=Exudation_Parameter(0);
[Mpar_H(1)]=Vegetation_Management_Parameter;
Mpar_H(1).Date_log = datenum(1923,1,15,12,0,0) ; %%%
Mpar_H(1).fract_log = 1; %%%
Mpar_H(1).fract_resprout = 0.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L = [NaN]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L= [NaN]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
%PLNR_L = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_L = [NaN];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_L= [NaN]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_L= [2.0 10.0];
aSE_L= [NaN]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_L= [NaN];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L =  [NaN]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_L = [NaN]; %% [°C] Cold Leaf SLed
drn_L=  [NaN]; %% turnover root  [1/d]
dsn_L= [NaN]; % normal transfer rate sapwood [1/d]
age_cr_L= [NaN]; %% [day] Critical Leaf Age
Bfac_lo_L= [NaN]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_L = [NaN]; %% Mean Temperature for Leaf onset
Tls_L = [NaN]; %% Mean Temperature for Leaf Shed
PAR_th_L= [NaN];
dmg_L= [NaN]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_L = [NaN];
Trr_L = [NaN]; %% Translocation rate [gC /m^2 d]
mjDay_L = [NaN]; %% Maximum Julian day for leaf onset
LDay_min_L =[NaN]; %% Minimum Day duration for leaf onset
LtR_L = [NaN]; %%% Leaf to Root ratio maximum
Mf_L= [NaN]; %% fruit maturation turnover [1/d]
Wm_L= [NaN] ; % wood turnover coefficient [1/d]
eps_ac_L = [NaN]; %% Allocation to reserve parameter [0-1]
Klf_L =[NaN]; %% Dead Leaves fall turnover [1/d]
fab_L = NaN; %% fraction above-ground sapwood and reserve
fbe_L =NaN; %% fraction below-ground sapwood and reserve
LDay_cr_L = [NaN]; %%%  Threshold for senescence day light [h]
ff_r_L= [NaN];
[ParEx_L(1)]=Exudation_Parameter(0);
[Mpar_L(1)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%% PRODUCTIVITY
%%% Maximum Rubisco Capacity
%[Vmax_H]=Maximum_Rubisco_Capacity(Sl_H,PLNR_H,Nl_H);
%[Vmax_L]=Maximum_Rubisco_Capacity(Sl_L,PLNR_L,Nl_L);
Vmax_H = [65]; %55
Vmax_L = [0]; %35
%[Amax_H]= MAX_PHOTOSYNTESIS(Vmax_H,Ca,CT_H,Tup_H,Tlow_H,FI_H,Oa); %% 17
%[Amax_L]= MAX_PHOTOSYNTESIS(Vmax_L,Ca,CT_L,Tup_L,Tlow_L,FI_L,Oa); %% 13
%%%%%%%%%%%%%%%%%%%
%Amax_L = [20]; %% 26
%Amax_H = [0]; %% 15
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%% Initial Condtion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%%%%%%%%%%%%%%%
%%% B1 Leaves - Grass  %%% B2 Sapwood  %%% B3 Fine Root  %%% B4 Carbohydrate Reserve
%%% B5 Fruit and Flower
%%% B6 Heartwood - Dead Sapwood
%%% B7 Leaves - Grass -- Standing Dead
LAI_L(1,:)=[0.0]; %
B_L(1,:,:)= [0 0 0 0 0 0 0 0]; %% 95 58 120 18
Rrootl_L(1,:)= [0] ;
PHE_S_L(1,:)=[0];
dflo_L(1,:)=[0];
AgeL_L(1,:)=[0];
e_rel_L(1,:)=[1];
hc_L(1,:) =[0.0]; %% 0.7
SAI_L(1,:) = [0.0]; %% 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAI_H(1,:)=[0];
%B_H(1,:,:)= [ 0 780 545 510 10 8073 18 0];
%B_H(1,:,:)= [ 0 780 545 510 10 12100 18 0];
B_H(1,:,:)= [ 0 780 545 510 10 0 18 0];
Rrootl_H(1,:)= [4000] ;
PHE_S_H(1,:)=[3];
dflo_H(1,:)=[0];
AgeL_H(1,:)=[0];
e_rel_H(1,:)=[1];
hc_H(1,:) =[22.0];
SAI_H(1,:) = [0.2];
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [400];
Preserve_H(1,:)= [50];
Kreserve_H(1,:)= [200];
FNC_H(1,:)=1;
NupI_H(1,:)= [0 0 0];
Nreserve_L(1,:)= [100];
Preserve_L(1,:)= [100];
Kreserve_L(1,:)= [100];
FNC_L(1,:)=1;
NupI_L(1,:)= [0 0 0];
RexmyI(1,:)= [0 0 0];
%%%
%%%%
if OPT_SoilBiogeochemistry == 1
    %%%%%
    %Wm_H= [1/22400]; %% 63 yr //  ; % wood turnover coefficient [1/d]
    Wm_H= [1/18624];%% 50.9 yr [1/17546]; %%% 48.04 yr
    B_H(1,1,6)= 8073; %% %%% 12552 with  48.04 yr    10657 with 50.9    with 63 yr  16458
    %%%%
    %Nreserve_H(1,:)= [1];
    %Preserve_H(1,:)= [1];
    %Kreserve_H(1,:)= [1];
    FNC_H(1,:)=[1 ];
    NupI_H(1,:)= [0.0232   0.003215   0.016318];
    RexmyI(1,:)= [0.082033252106261   0.086213955755071                   0];
    Nreserve_L(1,:)= [0]; Preserve_L(1,:)= [0]; Kreserve_L(1,:)= [0]; FNC_L(1,:)=[1];
    NavlI(1,:)=[0.298441431160444   0.183783586479934   0.338856740814781];
    %%%%
    %%% after spin-up
    B_H(1,1,:) = [0 555 447 396 6.5 6828 0.31 0];
    Nreserve_H(1,:)= 1.18; %
    Preserve_H(1,:)= 0.42;%
    Kreserve_H(1,:)= 6.94; %
    %%%%%%%
    %%%%% pre spin-up
    %B_H(1,1,:) = [0 596 463 425 6.9 7078 0.33 0];
    %Nreserve_H(1,:)= [1];
    %Preserve_H(1,:)= [1];
    %Kreserve_H(1,:)= [1];
    
    %%%% Spin-up
    
    P(1,:)=  1000*[0.110846823735348   0.142573118584581   0.023210840491207   0.810098548175676   0.270032849391892   0.031065897859492,...
        0.103953488264114   0.022376327271604   0.939310572365553   0.539620778986031   3.171787432562070   0.016418076171686,...
        0.012373580903029   0.000066551728836   0.000039663148618   0.000035797574885   0.000066105247696   0.023415088497764,...
        0.064073547275476   0.005970959380323   0.033835436488496   0.000035701647490   0.002906738108265   0.005220635088243,...
        0.001092159958430   0.350089306806673   0.004507012071765   0.009865049363892   0.000331994979927   0.001881329369643,...
        0.000120820409961   0.000013582874948   0.000010112740684   0.000003570164750   0.000658588672469   0.000372902506303,...
        0.000148554626951   0.075714069207306   0.001462225585361   0.001604226345466   0.000049799246989   0.000282199405446,...
        0.000188329845695   0.150000000000000   0.008275997578424   0.015029873994217   0.000000140884099   0.004263296130881,...
        0.002610317544122   0.000253661434585   0.009263534644226   0.000363353645055   0.000359565564357   0.009286353882914,...
        0.503286596584109];
    
    %%%%% Starting pre-spin-up
    % P(1,:)=  1000*[  0.115675907295798   0.144349510494787   0.023526254233599   0.837207887865445   0.279069295955147   0.032408354083219,...
    %0.106537243155629   0.023254797000246   1.086494181576372   0.553939149919292   3.240980016398846   0.016504547157517,...
    %0.012407712425145   0.000067145850001   0.000040108060153   0.000035599094042   0.000066846766921   0.024043405342063,...
    %0.066453270726474   0.005730762383795   0.032474320174841   0.000036368632713   0.003757460967849   0.005395339721800,...
    %0.001289826214977   0.396199144788530   0.004630618758761   0.010233455133466   0.000318639639192   0.001805647696027,...
    %0.000214275538558   0.000024709250181   0.000011966148106   0.000003636863272   0.000655276389674   0.000385381408700,...
    %0.000148464453988   0.079559849903521   0.001501479865338   0.001663842058318   0.000047795945879   0.000270847154404,...
    %0.000140938781628   0.150000000000000   0.005755259254664   0.015014058625103   0.000000147298169   0.004220472946484,...
    %0.002697669860900   0.000252816865016   0.009495683035372   0.000305898350106   0.000301231284293   0.008401304890628,...
    %0.503264878375030];
    
    
    %DepN=  0.0022; %   gN/m2 day
    %DepP= 2.3771e-006; %  gP/m2 day
    %DepK= 3.6986e-005; %%  gK/m2 day
    %FertN=0*ones(1,366);
    %FertP=0*ones(1,366);
    %FertK=0*ones(1,366);
    Upl=0.01; %% Soil production [mm/yr]
    HIST=0;
    %[B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl,HIST,FertN,FertP,FertK,DepN,DepP,DepK);
    [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl,HIST);
    B_IO.SC_par=[1 1 1 1];
    PHs = 4.1; %% 4.8
end
%%%
TBio_L=0;  %%[ton DM / ha ]
TBio_H=200; %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=10;  %% [mm/ m2 PFT];
Vl_H=10;  %% [mm/ m2 PFT];
Vx_L=0;   %% [mm/ m2 PFT];
Vl_L=0;   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%%%% Initialization OTHER
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+2;
Ts_under(1)=Ts(1);
Tdamp(1)=4.5;
Tdp(1,:)= 4.5*ones(1,ms);
TdpI_H(1,:)=4.5;
TdpI_L(1,:)=4.5;
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
O(1,:)= Ohy+0.15;
%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(cur_dir)
