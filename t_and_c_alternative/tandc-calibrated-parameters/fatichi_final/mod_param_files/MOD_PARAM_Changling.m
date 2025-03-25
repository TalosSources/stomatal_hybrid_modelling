%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETER FILE TEMPLATE FOR T&C %%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changling Grassland 

%Leymus chinensis, Phragmitis communis.
%Temperate, semi-arid continental monsoon climate
%Alkali-saline soil
%Heigth 0.7–0.8 m
%AG-Biomass  336.3 gC m-2 
% LAI 3.1
%%% GPP 415-570 gC m-2 yr-1
% ET 300-390 mm yr-1 
%Vegetation cover 70–80%
% 130–165 day frost free 
%pH >9.0 
%Soil Org. Mat. 1.5% -  35% clay, 45% silt, and 20% sand
%Tower height 6m 
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
zatm =6; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%
%%%% LAND COVER and VEGETATION COMPOSITION 
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare = 0.25; Ccrown = [0.74 0.01];      
%%% Rainfall disaggrgation information 
a_dis = NaN ;
pow_dis = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SOIL PARAMETERS 
Pcla= 0.35;     %
Psan= 0.20;      %              
Porg= 0.015;           
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
%%%%%%%%%%%%%
Zs= [0    10    20    50   100   200   300   400  500  600   800  1000  1200  1600  2000 ]; %% ms+1Zdes = 10;
Zdes = 10;               %  Soil depth [mm] relevant for evaporation (first layer)
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
WatFreez_Th = -8; %% [] Threshold for freezing water
dz_ice = 0.45; %% [mm/h] Freezing Layer depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%  ROOT PARAMETER 
ExEM = 1; %% Fraction of ectomychorrizae per area 
CASE_ROOT= 1;  %%% Type of Root Profile
ZR95_H = [0  0]; %% [mm]
ZR95_L = [700  700]; %% [mm]    [500 500]
ZR50_H = [NaN  NaN];
ZR50_L = [NaN  NaN];
ZRmax_H = [NaN  NaN];
ZRmax_L = [NaN  NaN];
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
Sp_LAI_H_In= [0.2         0.2]; %%[mm/LAI]
Sp_LAI_L_In= [0.2         0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [2  2]; %%[cm]
d_leaf_L= [0.5     0.5];  %% [cm]                            [0.7     1]
%% Veg Biochemical parameter
KnitH=[NaN  NaN]; %%% Canopy Nitrogen Decay
KnitL=[0.25 0.30];                  
%%%%%
mSl_H = [NaN  NaN];%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = [0  0]; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%
%%------
%%------
FI_H=[NaN  NaN];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[NaN  NaN]; %%[Pa] 
a1_H=[NaN  NaN];  %%% [-] WUE parameter 
go_H=[NaN  NaN];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[NaN  NaN]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[NaN  NaN];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[NaN  NaN]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[NaN  NaN]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[NaN  NaN]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------
FI_L=[0.081       0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]   [0.081       0.081]
Do_L=[800  800]; %%[Pa]     1000  1000]
a1_L=[7 7];  %%% [-] WUE parameter       [4 4]
go_L=[0.01        0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_L=[3  3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[0.649       0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[72  72]; %% [kJ / mol K]  entropy factor - Plant Dependent           
gmes_L=[Inf Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L=[1.9     2]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]    
%%%%%%%
Vmax_H = [0  0]; % [umol CO2 /m2 s] - Maximum Rubisco Capacity 
Vmax_L = [45  50]; % [umol CO2 /m2 s]                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hydraulic Parameters
Psi_sto_00_H = [NaN  NaN]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H = [NaN  NaN] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H = [NaN  NaN]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H = [NaN  NaN] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H = [NaN  NaN] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [NaN  NaN];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_H = [NaN  NaN] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [NaN  NaN];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [NaN  NaN]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [NaN  NaN]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata
Psi_sto_00_L= [-0.7        -0.7]; %% [MPa]  Water Potential at 2% loss conductivity             
Psi_sto_50_L = [-2.5        -5.5] ;%% [MPa]  Water Potential at 50% loss conductivity 
%%% Leaf
PsiL00_L = [-1.5        -1.4]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_L = [-3 -7] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [5  8] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [1200  1200];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_L = [6  6] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [80000  80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [-10  -9]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [150  150]; %%% [kg / m^3 sapwood MPa]
%% Growth Parameters
PsiG50_H= [NaN  NaN];  %%[MPa]
PsiG99_H= [NaN  NaN];  %%[MPa]
gcoef_H = [NaN  NaN]; % [gC/m2 day]
%%------
PsiG50_L= [-2.5        -1.4];        %%   
PsiG99_L= [-3        -5.5];          %%    
gcoef_L = [3.5         3.5]; % [gC/m2 day]
%%%
%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(0);
[PFT_opt_H(2)]=Veg_Optical_Parameter(0);
[PFT_opt_L(1)]=Veg_Optical_Parameter(13);
[PFT_opt_L(2)]=Veg_Optical_Parameter(9);
OM_H=[NaN  NaN];
OM_L=[1  1];
%%%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  VEGETATION PART %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% High Vegetation 
aSE_H=  [NaN  NaN]; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_H = [NaN  NaN]; %  [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_H = [NaN  NaN]; %[kgC/kgN ] Leaf Nitrogen Concentration
r_H =  [NaN  NaN];  %% respiration rate at 10?[gC/gN d ]  [0.066 -0.011]
gR_H=  [NaN  NaN];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
dd_max_H= 1./[NaN  NaN];%%%  [1/d]  death maximum for drought
dc_C_H =  1./[NaN  NaN]; %%[Factor of increasing mortality for cold]
Tcold_H = [NaN  NaN]; %% [] Cold Leaf Shed
drn_H=  1./[NaN  NaN]; %% turnover root  [1/d]
dsn_H= 1./[NaN  NaN]; % normal transfer rate sapwood [1/d]
age_cr_H= [NaN  NaN]; %% [day] Critical Leaf Age
Trr_H = [NaN  NaN]; %% Translocation rate [gC /m^2 d]
LtR_H = [NaN  NaN]; %%% Leaf to Root ratio maximum
Mf_H= 1./[NaN  NaN]; %% fruit maturation turnover [1/d]
Wm_H= [NaN  NaN] ; % wood turnover coefficient [1/d]
eps_ac_H = [NaN  NaN]; %% Allocation to reserve parameter [0-1]
Klf_H = 1./[NaN  NaN]; %% Dead Leaves fall turnover [1/d]
fab_H = [NaN  NaN]; %% fraction above-ground sapwood and reserve
fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve
ff_r_H= [NaN  NaN]; %% Reference allocation to Fruit and reproduction
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
[Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
%%%% Phenology 
Bfac_lo_H= [NaN  NaN]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN  NaN] ;  % Not-used 
Tlo_H = [NaN  NaN]; %% Mean Temperature for Leaf onset
Tls_H = [NaN  NaN]; %%% Not-used 
PAR_th_H= [NaN  NaN]; %% Light Phenology Threshold 
dmg_H= [NaN  NaN]; %%%  Day of Max Growth
LAI_min_H = [NaN  NaN];
mjDay_H = [NaN  NaN]; %% Maximum Julian day for leaf onset
LDay_min_H =[NaN  NaN]; %% Minimum Day duration for leaf onset
LDay_cr_H = [NaN  NaN]; %%%  Threshold for senescence day light [h]
[ParEx_H(1)]=Exudation_Parameter(0);
[Mpar_H(1)]=Vegetation_Management_Parameter;
[ParEx_H(2)]=Exudation_Parameter(0);
[Mpar_H(2)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Low Vegetation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
aSE_L=  [2  0]; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_L = [0.018   0.010]; %  [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005     [0.02   0.008]    
Nl_L = [35  45]; %[kgC/kgN ] Leaf Nitrogen Concentration
r_L =  [0.045  0.056];  %% respiration rate at 10?[gC/gN d ]  [0.066 -0.011]
gR_L=  [0.25   0.25];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]               
dd_max_L= 1./[25  365];%%%  [1/d]  death maximum for drought
dc_C_L =  1./[52  182]; %%[Factor of increasing mortality for cold]
Tcold_L = [6  -35]; %% [][12  -25]     Cold Leaf Shed                                         
drn_L=  1./[1200   900]; %% turnover root  [1/d]
dsn_L= 1./[365  600]; % normal transfer rate sapwood [1/d]
age_cr_L= [200  730]; %% [day] Critical Leaf Age         [200  365]       
Trr_L = [1.6   1.2]; %% Translocation rate [gC /m^2 d]    [1.0   0.5]
LtR_L = [0.7    0.5]; %%% Leaf to Root ratio maximum
Mf_L= 1./[30  30]; %% fruit maturation turnover [1/d]
Wm_L= [0  0] ; % wood turnover coefficient [1/d]
eps_ac_L = [1     0.8]; %% Allocation to reserve parameter [0-1]
Klf_L = 1./[25  25]; %% Dead Leaves fall turnover [1/d]
fab_L = [0        0.5]; %% fraction above-ground sapwood and reserve
fbe_L = 1-fab_L; %% fraction below-ground sapwood and reserve
ff_r_L= [0.1         0.1]; %% Reference allocation to Fruit and reproduction
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
[Stoich_L(2)]=Veg_Stoichiometric_Parameter(Nl_L(2));
%%%% Phenology 
Bfac_lo_L= [0.96 0.96]; %% Leaf Onset Water Stress  
Bfac_ls_L= [NaN  NaN] ; % Not-used 
Tlo_L = [3.8  2.5]; %%     [-4  0]  Mean Temperature for Leaf onset    
Tls_L = [NaN  NaN]; %% Not-used 
PAR_th_L= [Inf Inf]; %% Light Phenology Threshold 
dmg_L= [30  20]; %%%  Day of Max Growth                     
LAI_min_L = [0.05       0.001];
mjDay_L = [366  366]; %% Maximum Julian day for leaf onset
LDay_min_L =[13.2  13.2]; %% Minimum Day duration for leaf onset   
LDay_cr_L = [12.7  12.7]; %%%  Threshold for senescence day light [h]    
[ParEx_L(1)]=Exudation_Parameter(0);
[Mpar_L(1)]=Vegetation_Management_Parameter;
[ParEx_L(2)]=Exudation_Parameter(0);
[Mpar_L(2)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
LAI_H(1,:)=[0  0]; %
B_H(1,:,:)= [0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0]; %% 
Rrootl_H(1,:)= [0  0] ;
PHE_S_H(1,:)=[0  0];
dflo_H(1,:)=[0  0];
AgeL_H(1,:)=[0  0];
e_rel_H(1,:)=[1  1];
hc_H(1,:) =[0  0]; %%
SAI_H(1,:) = [0  0]; %% 
%%%%%%%%%%%%%%%%%%%%%
LAI_L(1,:)=[ 0    0];                             
B_L(1,:,:)= [0  0 255 201 0 0  118 0; 
             290 324 580 253 0 0 0 0];       
Rrootl_L(1,:) = [2667 6073] ;      
PHE_S_L(1,:)=[1 1];
dflo_L(1,:)=[0 0];
AgeL_L(1,:)=[0  723];              
e_rel_L(1,:)=[1 1];                          
hc_L(1,:) =[0.02   0.1];         
SAI_L(1,:) = [0.001 0.05];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [0  0];
Preserve_H(1,:)= [0  0];
Kreserve_H(1,:)= [0  0];
FNC_H(1,:)=[1  1];
NupI_H(1,:,:)= [0 0 0 ;0 0 0];
Nreserve_L(1,:)= [1000   1000 ];
Preserve_L(1,:)= [1000    1000 ];
Kreserve_L(1,:)= [1000   1000 ];
FNC_L(1,:)=[1  1];
NupI_L(1,:,:)= [0 0 0 ;0 0 0];
RexmyI(1,:)= [0 0 0];
%%%
%%%%
if OPT_SoilBiogeochemistry == 1
    %%%%%
    %%%%
end
%%%
TBio_L=[1  1];  %%[ton DM / ha ]
TBio_H=[1  1]; %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=[1  1];  %% [mm/ m2 PFT];
Vl_H=[1  1];  %% [mm/ m2 PFT];
Vx_L=[1  1];   %% [mm/ m2 PFT];
Vl_L=[1  1];   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%%%% Initialization OTHER
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+1;
Tdamp(1)=-3;
Tdp(1,:)= [  -18.6341  -18.2454  -17.4715  -15.9921  -13.6297  -11.8072  -10.6781   -9.7444   -8.8392   -7.4452   -5.5993   -3.8837    -1.5142   -0.3221]; 
TdpI_H(1,:)= 0;
TdpI_L(1,:)= -5;  
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
O(1,:)=[   0.2638    0.2638    0.2638    0.2638    0.2638    0.2638    0.2638    0.2638    0.2638    0.2660    0.2727    0.2818 0.2076    0.2580];
Oice(1,:)= [  0.0495    0.0489    0.0473    0.0389    0.0291    0.0480    0.0650    0.0810    0.1236    0.0876    0.0929    0.0941  0.0563    0.0093];  
%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
Vice(1,:) = Oice(1,:).*dz;
%%%%%%%%%%%%%%%%%
cd(cur_dir)
%%%%%%%%%%%%%%%%%