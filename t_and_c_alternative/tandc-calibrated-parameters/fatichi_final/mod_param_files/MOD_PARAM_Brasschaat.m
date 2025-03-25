

%Tower height 40m 
%This forest consists of 117
%patches of a wide variety of species combinations with
%approximately 80% of the surface area of the forest
%covered by trees, and 16% comprised of grasses
%(predominately Molinia caerulea (L.). Although several tree species
%are present, pedunculate oak (Quercus robur L.; 80% of
%the deciduous species, mostly planted in 1936) and
%Scots pine (Pinus sylvestris L.; 80% of the coniferous
%species, mostly planted in 1929) are the dominant overstory species
%mixed coniferous/deciduous forest. All undergrowth was completely removed in 1993. 

%GPP 1150 g C m/2 per year. 



%%%% Deciduous Broadleaf Forest
%%%Tower height 39m
%Scots pine forest, planted in 1929
%Scots pine (Pinus sylvestris L., 80% of the coniferous species) and pedunculate oak (Quercus robur L., 80% of the deciduous species) dominate the canopy composition. The undergrowth is dominated by black cherry (Prunus serotina Ehrh.), rhododendron (Rhododendron ponticum L.) and grass (Molinia caerulea L. Moench), and an extensive moss cover (dominated by Hypnum cupressiforme Hedw.)
%%%// 30% evergreen 30% deciduous 25% grass 15% bare soil 
%Near urban area on sandy soil  sand deposits covering a 0.5 m thick loamy  clay layer at  1.75 m of depth
%Rain 750
%Mean annual temp 9.8°C
%%% Trre height 21m 
%%% Root depth 160 cm 
%%% Phen 109-117 begin  291-318 end 
%%% GPP 1100 -1200 
%LAI	3  (1.7 - 3.6 winter summer) 
%A groundwater depth usually is between 1.2 and 1.5 m

%%%%%%%%%%%% VEGETATION COMPOSITION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOIL AND HYDROLOGICAL PARAMETER
%%%%%%%%%%%%%%%%%
cur_dir=cd;
cd(Directory)
%%%%
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
Ared = 1; %% 1-Frock, Reduction in Area due to rock content 
aR =1; %%% anisotropy ratio
%Kh=Ks*aR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsize= 1; %%[m^2];
aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
Kbot = NaN; %% [mm/h] Conductivity at the bedrock layer
Krock = NaN; %% [mm/h] Conductivity of Fractured Rock 
zatm = 39; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%%%%%% VEG. SPECIES  --- Low Grasses -- High Decidous
%%%% LAND COVER PARTITION
Cwat = 0; Curb = 0 ; Crock = 0;
Cbare = 0; Ccrown = [0.60 0.40];
%%%%%%%%%% SOIL INPUT - 
Pcla= 0.15; 
Psan= 0.80;
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
nVG=L+1;
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
Ks_Zs(14:end)=Ks_Zs(14:end)/100; 
%%%%%%%%%%%%%
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
clear Oss_Lp Owp_Lp 
Oice = 0;
%%%%%%%%%%%%%
Zs= [0 10 20 50 100 200 300 400 600 800 1000 1200 1500 1750 2000 2250 2500 ]; %% ms+1
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
ExEM = 0.5;
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS VEGETATION
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile 
ZR95_H = [1600 0]; %% [mm]
ZR95_L = [0 1600]; %% [mm]
ZR50_H = [NaN NaN];
ZR50_L = [NaN NaN]; 
ZRmax_H = [NaN NaN]; 
ZRmax_L = [NaN NaN]; 
%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
%5 Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_L_In= [0.2 0.2]; %%[mm/LAI]
Sp_LAI_H_In= [0.1 0.1]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [0.25 0.25]; %%[cm]
d_leaf_L= [3.5 3.5];  %% [cm]
%%%%%%%% Biochemical parameter
%KnitH=0.35; %%% Canopy Nitrogen Decay
KnitH=[0.35 0.35]; %%% Canopy Nitrogen Decay
KnitL=[0.30 0.30] ; 
mSl_H = [0.0 0.0];%% [m2 PFT /gC]  Linear increase in Sla with LAI 
mSl_L = [0.0 0.0]; 
%%------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_H=[0.081 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[700 700]; %%[Pa]
a1_H=[6 6];
go_H=[0.01 0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3 3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance 
rjv_H=[2.0 2.0]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
FI_L=[0.081 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[1000 1000]; %%[Pa]
a1_L=[7 7];
go_L=[0.01 0.01];% % [mol / s m^2] minimum Stomatal Conductance
CT_L=[3 3];  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[0.649 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[Inf Inf];
rjv_L= [2.0 2.0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Matric Potential
%Pss_H = [1000 1000]; 
%Pwp_H = [3200 3200]; %%% [kPa]
%Pss_L = [800 800]; 
%Pwp_L = [2800 2800]; %%% [kPa]
%%% Stomata
Psi_sto_00_H =  [-0.8 -0.8]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H =  [-2.0 -2.0] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H =  [-1.0 -1.0] ;%%[MPa]  Water Potential at 50% loss conductivity
PsiL50_H =  [-2.6 -2.6]; %% [MPa]  Water Potential at 2% loss conductivity
Kleaf_max_H = [ 20 20] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200 1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = [6.0 6.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-6.5 -6.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [80 80]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata %% -1.4 4.0 
Psi_sto_00_L =  [-0.8 -0.8]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_L =  [-2.2 -2.2] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L =  [-1.2 -1.2]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_L =  [-2.5 -2.5] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [ 8 8] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [1200 1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L = [6.0 6.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [-4.5 -6.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [150 150]; %%% [kg / m^3 sapwood MPa]

%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L); 
%%%% Growth Parameters 
PsiG50_H= [-0.8;-0.8];  %%[MPa]
PsiG99_H= [-2.0; -2.0] ; %%[MPa]
gcoef_H = [3.5 3.5]; % [gC/m2 day]
%%------  
PsiG50_L= [-0.8;-0.8];  %%[MPa]
PsiG99_L= [-2.2; -2.2] ; %%[MPa]
gcoef_L = [3.5 3.5]; % [gC/m2 day]

%%%%%%%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(1);%could be 1
[PFT_opt_L(1)]=Veg_Optical_Parameter(7);
[PFT_opt_H(2)]=Veg_Optical_Parameter(1);%could be 1
[PFT_opt_L(2)]=Veg_Optical_Parameter(7);
OM_H=[1 1];
OM_L=[1 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.011 0]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
%Sl_H = [0.0112 0];
Nl_H= [42 0]; %[gC/gN ] Leaf Carbon-Nitrogen ratio 
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
[Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
r_H = [0.050 0];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_H= [0.25 0]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_H= [2.0 10.0];
aSE_H= [0 0]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops 
dd_max_H= [1/200 0]; %%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  [78/365 0]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H = [-20.0 0]; %% [°C] Cold Leaf Shed
drn_H=  [1/700 0]; %% turnover root  [1/d]
dsn_H= [1/500 0]; % normal transfer rate sapwood [1/d] 
age_cr_H= [650 0]; %% [day] Critical Leaf Age
Bfac_lo_H= [0.99 0]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [6.5 0]; %% Mean Temperature for Leaf onset
Tls_H = [NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_H= [NaN NaN]; 
dmg_H= [35 0]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.001 0];
Trr_H = [0.5 0]; %% Translocation rate [gC /m^2 d]
mjDay_H = [250 0]; %% Maximum Julian day for leaf onset
LDay_min_H =[13.1 0]; %% Minimum Day duration for leaf onset
LtR_H = [0.8 0]; %%% Leaf to Root ratio maximum 
Mf_H= [1/50 0]; %% fruit maturation turnover [1/d]
Wm_H= [0 0] ; % wood turnover coefficient [1/d]
eps_ac_H = [0.3 0]; %% Allocation to reserve parameter [0-1] 
LDay_cr_H = [10.8 0]; %%%  Threshold for senescence day light [h]
Klf_H =[1/25 0]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74 0]; %% fraction above-ground sapwood and reserve
fbe_H = [0.26 0]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1 0.1]; %% Reference allocation to Fruit and reproduction 
[ParEx_H(1)]=Exudation_Parameter(1); 
[ParEx_H(2)]=Exudation_Parameter(1); 
[Mpar_H(1)]=Vegetation_Management_Parameter;
[Mpar_H(2)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L = [0 0.028]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
%Sl_L = [0 0.0319]; 
Nl_L= [0 30]; %[kgC/kgN ] Leaf Nitrogen Concentration
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
[Stoich_L(2)]=Veg_Stoichiometric_Parameter(Nl_L(2));
r_L =  [0 0.050];  %% [0.066 -0.011] respiration rate at 10° [gC/gN d ]
gR_L=  [0 0.25];% [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_L=  [ NaN NaN];
aSE_L=  [1 1]; %%%%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_L= [0 1/200];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L =  [0 55/365]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_L = [0 5.0]; %% [°C] Cold Leaf Shed
drn_L=  [0 1/600]; %% turnover root  [1/d]
dsn_L= [0 1/630]; % normal transfer rate sapwood [1/d] 
age_cr_L= [0 140]; %% [day] Critical Leaf Age
Bfac_lo_L= [0 0.99]; %% Leaf Onset Water Stress
Bfac_ls_L= [0 NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_L = [0 7.2]; %% Mean Temperature for Leaf onset 11.3 
Tls_L = [NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_L= [NaN NaN]; 
dmg_L= [0 28]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_L = [0 0.05];
Trr_L = [0 7.0]; %% Translocation rate [gC /m^2 d]
mjDay_L = [0 250]; %% Maximum Julian day for leaf onset
LDay_min_L =[0 13.3]; %% Minimum Day duration for leaf onset 
LtR_L = [0 0.7]; %%% Leaf to Root ratio maximum 
Mf_L= [0 1/50]; %% fruit maturation turnover [1/d]
Wm_L= [0 0] ; % wood turnover coefficient [1/d]
eps_ac_L = [0 1]; %% Allocation to reserve parameter [0-1] 
LDay_cr_L = [0 10.8]; %%%  Threshold for senescence day light [h]
Klf_L =[0 1/28]; %% Dead Leaves fall turnover [1/d]
fab_L = [0 0.74]; %% fraction above-ground sapwood and reserve
fbe_L = [0 0.26]; %% fraction below-ground sapwood and reserve
ff_r_L= [0.1 0.1];
[ParEx_L(1)]=Exudation_Parameter(0); 
[ParEx_L(2)]=Exudation_Parameter(0); 
[Mpar_L(1)]=Vegetation_Management_Parameter;
[Mpar_L(2)]=Vegetation_Management_Parameter;
%%%%%%%%%%%% PRODUCTIVITY
%%% Maximum Rubisco Capacity
%[Vmax_H]=Maximum_Rubisco_Capacity(Sl_H,PLNR_H,Nl_H);
%[Vmax_L]=Maximum_Rubisco_Capacity(Sl_L,PLNR_L,Nl_L);
%Vmax_H = [80]; %55
Vmax_H = [38 0]; %55
Vmax_L = [0 48];%[45]; %35
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
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
%%%%%%%
LAI_H(1,:)=[2.7 0]; %
B_H(1,:,:)= [271 385 339 279 2 0 10 0; 
    0 0 0 0 0 0 0 0 ]; %%
Rrootl_H(1,:)= [0 0] ;
PHE_S_H(1,:)=[4 0];
dflo_H(1,:)=[0 0];
AgeL_H(1,:)=[663 0];
e_rel_H(1,:)=[0 0];
hc_H(1,:) =[21 0]; %%
SAI_H(1,:) = [0.02 0]; %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LAI_L(1,:)=[0 0.0]; %
B_L(1,:,:)= [0 0 0 0 0 0 0 0
    0 324 171 234 2 0 13 0 ]; %% 
Rrootl_L(1,:) = [1000 1000] ; 
PHE_S_L(1,:)=[0 4];
dflo_L(1,:)=[0 0];
AgeL_L(1,:)=[0 0];
e_rel_L(1,:)=[0 0];
hc_L(1,:) =[0 10]; %%
SAI_L(1,:) = [0 0.01]; %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass 
%%%%%%%%%%
Nreserve_H(1,:)= [1000 1000 ];
Preserve_H(1,:)= [100 100];
Kreserve_H(1,:)= [100 100];
FNC_H(1,:)=[ 1 1 ];
NupI_H(1,:,:)= [0 0 0 ; 0 0 0];
Nreserve_L(1,:)= [1000 1000];
Preserve_L(1,:)= [100 100];
Kreserve_L(1,:)= [100 100];
FNC_L(1,:)=[1 1];
NupI_L(1,:,:)= [0 0 0; 0 0 0];
RexmyI(1,:)= [0 0 0];
%%%%%%%%
%%%%%%%
if OPT_SoilBiogeochemistry == 1
    %%%%%
end
%%%%%%%%
TBio_L=[50 50];  %%[ton DM / ha ]
TBio_H=[200 200]; %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=[0 0];  %% [mm/ m2 PFT];
Vl_H=[100 100];  %% [mm/ m2 PFT];
Vx_L=[0 0];   %% [mm/ m2 PFT];
Vl_L=[100 100];   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+2;
Tdamp(1)=Ta(1); 
Tdp(1,:)= Ta(1)*ones(1,ms);
Sdp(1,:)= (0.001*dz).*cv_s.*(Tdp(1,:)); %% [J °C/m^2 K]
TdpI_H(1,:)=5;
TdpI_L(1,:)=5;
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
Ws_under(1,:)=1; 
%%%%%%%%%%%%%% Volume [mm]
O(1,:)= Ohy+0.25;
%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
cd(cur_dir)
%%%%%%%%%%%%%%%%%