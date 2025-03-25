%%%%%%%%%%%% VEGETATION COMPOSITION 

%%%%%%% GPP = 1600 gC/m2 year
%%% NPP = 700 gC/m2 
%%% LAI = 4.0 
%%% Stand Age = 20 years   5 kgC m-2 
%%% Fine root biomass 400 gDM m-2 

%%% Cambisol, sandy-clay Turkey oak forest
%%% The bedrock is of volcanic origin and the soil is a 90 cm deep Luvisol
%%% 40-50% Sand , 30-50% Clay 

%16–20 m high main tree canopy, largely made up of Q. cerris, which includes also
%Quercus pubescens Willd., Fraxinus ornus L., Ulmus minor Mill., Ostrya carpinifolia
%Scop., Acer monspessulanum L. and A. campestre L. An understorey layer typically
%composed of Phillyrea latifolia L., Crataegus spp., Ruscus aculeatus L., Cornus spp.,
%Prunus spinosa L., Sorbus spp. completes the forest vegetation

%Tower height 20 m 5 m above canopy 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOIL AND HYDROLOGICAL PARAMETER
%%%%%%%%%%%%%%%%%
cur_dir=cd;
cd(Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
zatm = 20; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%%% LAND COVER PARTITION
Cwat = 0; Curb = 0.0 ; Crock = 0.0;
Cbare = 0; Ccrown = [1.0];
%%%%%%%%%% SOIL INPUT %%%  sandy-clay
Pcla= 0.30; 
Psan= 0.45;
Porg= 0.03; 
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
%%% Measured %%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%%%%%%%%%%%%%%%
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
clear Oss_Lp Owp_Lp 
Oice = 0;
%%%%%%%%%%%%%%
Zs = [0 10 50 100 150 200 300 400 500 600 800 1000 1200 1500 2000]; %%% [ms+1]
Zdes = 10;
Zinf = 50; 
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
%%%%%%%%%%%%%%% OTHER PARAMETER
In_max_urb=5; In_max_rock=0.1; %% [mm]
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
ExEM = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS VEGETATION
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile 
ZR95_H = [1500]; %% [mm]
ZR95_L = [0]; %% [mm]
ZR50_H = NaN;
ZR50_L = NaN; 
ZRmax_H = NaN;
ZRmax_L = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_L_In= [0.1]; %%[mm/LAI]
Sp_LAI_H_In= [0.15]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [4]; %%[cm]
d_leaf_L= [0.8];  %% [cm]
%%%%%%%% Biochemical parameter
KnitH=0.25; %%% Canopy Nitrogen Decay
KnitL=0.2; 
mSl_H = 0.0;%% [m2 PFT /gC]  Linear increase in Sla with LAI 
mSl_L = 0.0; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree 
%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Photosynthesis Parameter
%%------
FI_H=[0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[600]; %%[Pa]
a1_H=[7];
go_H=[0.001];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649];  %% [kJ/mol] Activation Energy - Plant Dependent 
Ha_H =[76]; %% [kJ / mol K]  entropy factor - Plant Dependent 
gmes_H=[Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance 
rjv_H=[2.2]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FI_L=[0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[1400]; %%[Pa]
a1_L=[7];
go_L=[0.01];% % [mol / s m^2] minimum Stomatal Conductance
CT_L=[3];  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[0.656];  %% [kJ/mol] Activation Energy - Plant Dependent 
Ha_L =[55]; %% [kJ / mol K]  entropy factor - Plant Dependent 
gmes_L=[Inf];
rjv_L= [2.6];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pss_H = [1000]; Pwp_H = [3500]; %%% [kPa]
%Pss_L = [500]; Pwp_L = [3500]; %%% [kPa]
Psi_sto_00_H =  -0.5; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H =  -3.2 ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H =  -1.0; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H =  -3.8 ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H = 10 ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = 1200;  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = 15.0 ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = 80000;  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = -9.0; %%[MPa]  Water Potential at 50% loss conductivity
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
PsiG50_H= -0.5;  %%[MPa]
PsiG99_H= -3.2;  %%[MPa]
gcoef_H = 4.5; % [gC/m2 day]
%%------  
PsiG50_L= -1.45; 
PsiG99_L= -4.0; 
gcoef_L = 3.5; % [gC/m2 day]
%%%%%%%% Vegetation Optical Parameter 
[PFT_opt_H(1)]=Veg_Optical_Parameter(7);
[PFT_opt_L(1)]=Veg_Optical_Parameter(13);  
OM_H=[1 ];
OM_L=[1 ];
%%%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.020]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [28]; %[gC/gN ] Leaf Carbon-Nitrogen ratio 
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
r_H = [0.032];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_H= [0.25]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_H= [2.0 10.0];
aSE_H= [1]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops 
dd_max_H= [1/100];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  36/365; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H = [6.5]; %% [°C] Cold Leaf Shed
drn_H=  [1/1400]; %% turnover root  [1/d]
dsn_H= [1/1100]; % normal transfer rate sapwood [1/d] 
age_cr_H= [180]; %% [day] Critical Leaf Age
Bfac_lo_H= [0.99]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [12.5]; %% Mean Temperature for Leaf onset 11.3 
Tls_H = [NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_H= [NaN];
dmg_H= [30]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.05];
Trr_H = [7.0]; %% Translocation rate [gC /m^2 d]
mjDay_H = [180]; %% Maximum Julian day for leaf onset
LDay_min_H =[12.0]; %% Minimum Day duration for leaf onset 
LtR_H = [0.6]; %%% Leaf to Root ratio maximum 
Mf_H= [1/80]; %% fruit maturation turnover [1/d]
Wm_H= [0] ; % wood turnover coefficient [1/d]
eps_ac_H = [1]; %% Allocation to reserve parameter [0-1] 
LDay_cr_H = [11.8]; %%%  Threshold for senescence day light [h]
Klf_H =[1/28]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74]; %% fraction above-ground sapwood and reserve
fbe_H = [0.26]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1]; %% Reference allocation to Fruit and reproduction 
[ParEx_H(1)]=Exudation_Parameter(0);
[Mpar_H(1)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L = [0.035]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L= [23]; %[kgC/kgN ] Leaf Nitrogen Concentration
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
r_L =  [ 0.045];  %% [0.066 -0.011] respiration rate at 10° [gC/gN d ]
gR_L=  [ 0.25];% [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_L=  [ NaN NaN];
aSE_L=  [2]; %%%%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_L=  [1/250];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L =  [6/365]; %% [Factor of increasing mortality for cold]
Tcold_L =  [0]; %% [°C] Cold Leaf Shed
drn_L=   [1/450]; %% senescence rate leaf  [1/d]
dsn_L=  [1/365]; % [ normal transfer rate sapwood [1/d] ]
age_cr_L= [ 180]; %% [day] Critical Leaf Age
Bfac_lo_L=  [0.99]; %% PHENOLOGY Leaf Onset Water Stress [0-1]
Bfac_ls_L=  [ NaN]; %% PHENOLOGY Leaf Shed Water Stress [0-1]
Tlo_L =  [1.0]; %% Mean Temperature for Leaf onset
Tls_L =  [ NaN]; %% Mean Temperature for Leaf Shed
PAR_th_L= [NaN];
dmg_L=  [ 20]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_L =  [0.1];
Trr_L =  [ 2.0]; %% Translocation rate [gC /m^2 d]
mjDay_L = [ 250]; %% Maximum Julian day for leaf onset
LDay_min_L =[10.7]; %% Minimum Day duration for leaf onset
LtR_L = [0.35]; %%% Leaf to Root ratio maximum 
Mf_L= [1/50]; %% fruit maturation turnover [1/d]
Wm_L= [0] ; % wood turnover coefficient [1/d]
eps_ac_L = [0.2]; %% Allocation to reserve parameter [0-1] 
LDay_cr_L = [10.7]; %%%  Threshold for senescence day light [h]
Klf_L =[1/50]; %% Dead Leaves fall turnover [1/d]
fab_L = 0; %% fraction above-ground sapwood and reserve
fbe_L =1; %% fraction below-ground sapwood and reserve
ff_r_L= [NaN];
[ParEx_L(1)]=Exudation_Parameter(0); 
[Mpar_L(1)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%% PRODUCTIVITY
%%% Maximum Rubisco Capacity
%[Vmax_H]=Maximum_Rubisco_Capacity(Sl_H,PLNR_H,Nl_H);
%[Vmax_L]=Maximum_Rubisco_Capacity(Sl_L,PLNR_L,Nl_L);
Vmax_H = [72]; %55
Vmax_L = [0]; %35
%[Amax_H]= MAX_PHOTOSYNTESIS(Vmax_H,Ca,CT_H,Tup_H,Tlow_H,FI_H,Oa); %% 17
%[Amax_L]= MAX_PHOTOSYNTESIS(Vmax_L,Ca,CT_L,Tup_L,Tlow_L,FI_L,Oa); %% 13
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
 [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end 
Lmax_day = max(L_day); 
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%%%%%%%%%%%%%%%
LAI_H(1,:)=[0.8]; %
B_H(1,:,:)= [34 691 382 470 8 0 141 0]; %%
Rrootl_H(1,:)= [4000] ;
PHE_S_H(1,:)=[4];
dflo_H(1,:)=[238];
AgeL_H(1,:)=[206];
e_rel_H(1,:)=[1];
hc_H(1,:) =[16]; %% 0.7
SAI_H(1,:) = [0.2]; %% 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LAI_L(1,:)=[0.0];
B_L(1,:,:)= [17 40 40 25 0 0 0 0 ];
Rrootl_L(1,:)= [0] ;
PHE_S_L(1,:)=[1 ];
dflo_L(1,:)=[0  ];
AgeL_L(1,:)=[89];
e_rel_L(1,:)=[1];
hc_L(1,:) =[0.5];
SAI_L(1,:) = [0.001];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [1000];
Preserve_H(1,:)= [100];
Kreserve_H(1,:)= [100];
FNC_H(1,:)=1;
NupI_H(1,:,:)= [0 0 0];
Nreserve_L(1,:)= [1000];
Preserve_L(1,:)= [100];
Kreserve_L(1,:)= [100];
FNC_L(1,:)=1;
NupI_L(1,:,:)= [0 0 0];
RexmyI(1,:)= [0 0 0];
%%%%%%%%
%%%%%%
if OPT_SoilBiogeochemistry == 1

    %%%%% 
    %%%% 
    Wm_H= [1/24090]; %%   ; % wood turnover coefficient [1/d]
    %B_H(1,1,6)= 13500;
    %%%
    B_H(1,1,6)= 2600;

    
    Nreserve_H(1,:)= [2.0];%
    Preserve_H(1,:)= [0.5];%
    Kreserve_H(1,:)=  [2];%
    FNC_H(1,:)=[1 ];
    NupI_H(1,:)= [   0.016687982   0.00079064   0.0150349];
    RexmyI(1,:)= [   0.0264955   0.14992                  0];
    Nreserve_L(1,:)= [0]; Preserve_L(1,:)= [0]; Kreserve_L(1,:)= [0]; FNC_L(1,:)=[1];
    NavlI(1,:)=[  0.09843   0.00417   0.349108398];
    
%     %%%% Pre-spin up
%      P(1,:)=  1000*[  0.045108810926662   0.119656974110512   0.019844548675185   0.935890755113580   0.311963585037855   0.013874482957307,...
%    0.118454743986182   0.029984137628689   1.755039811109003   0.842797349105752   5.736030981748688   0.014060623683752,...
%    0.011248498946971   0.000079532039562   0.000054059253402   0.000033138349818   0.000090098755670   0.021360889652207,...
%    0.083597790867539   0.023261977278759  -0.000000000346416   0.001458050623250   0.001679387762311   0.006462102832926,...
%    0.001248655278782   0.609905557553590   0.004118163825334   0.012881226692595   0.001290616520720  -0.000000000170436,...
%    0.000333041934598   0.000030681211077   0.000012047009233   0.000143805062325   0.000119956268736   0.000461578773781,...
%    0.000089189662770   0.109453555722603   0.001333063914457   0.002093974481735   0.000193592478108  -0.000000000170436,...
%    0.000027863745208   0.150000000000000   0.000585220935043   0.015007688992453   0.000000122966693   0.000839693881155,...
%    0.003231051416463   0.000196397928094   0.007921776791966   0.000085765593851   0.000073760688664   0.001659612619917,...
%    0.503239744227819]; 

       %%%% Post

      P(1,:)=  1000*[0.029130475033903   0.061768797541223   0.010192747738331   0.272420896713590   0.090806965571197   0.017461218718079,...
   0.051724893493591   0.009963442485573   1.715730350953937   0.862301078187807   5.955043538278672   0.020173747548346,...
   0.014835564303744   0.000059091534248   0.000040957775127   0.000024621472603   0.000068262958545   0.014605458163721,...
   0.057417888747204   0.054775150660963                   0   0.001471194313562   0.001125116547542   0.001881001429689,...
   0.000796473130667   0.628575818891257   0.002812111299033   0.008845849402429   0.003046427492144                   0,...
   0.000067536566845   0.000006180965362   0.000009725500415   0.000145119431356   0.000057386864344   0.000134357244978,...
   0.000041401340087   0.132046589205019   0.000910641141427   0.001437992796506   0.000456964123822                   0,...
   0.000002408832062   0.150000000000000   0.000421170552423   0.015007786730114   0.000000100482636   0.001295220226694,...
   0.000940500714845   0.000158074147886   0.009520255406939   0.000088982012583   0.000080545891814   0.002219561490039,...
   0.503242803005020];
    
    Upl=0.01; %% Soil production [mm/yr]
    HIST=0;
    %[B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl,HIST,FertN,FertP,FertK,DepN,DepP,DepK);
    [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl,HIST);
    B_IO.SC_par=[1 1 1 1];
    PHs=5.7;  %

end
%%%
%%%%%%%%%%%%%%%%%
Vx_H=10;  %% [mm/ m2 PFT];
Vl_H=10;  %% [mm/ m2 PFT];
Vx_L=0;   %% [mm/ m2 PFT];
Vl_L=0;   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%%%%%%% Initial Conditions
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1);
Ts(1)=Ta(1)+2;
Ts_under(1)=12; 
Tdamp(1)=11; 
Tdp(1,:)= 11*ones(1,ms);
TdpI_H(1,:)=11;
TdpI_L(1,:)=11;
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
Ws_under(1)=1; 
%%%%%%%%%%%%%% Volume [mm]
%%%%%%%%%%%%%%%%%%%%%
%ZWT(1)= 1999; %% [mm]
O(1,:)= [ 0.3412    0.3476    0.3546    0.3572    0.3580    0.3584    0.3589    0.3598    0.3606    0.3613    0.3609    0.3600    0.3596    0.3612];
%dz(ms)=ZWT(1)-Zs(ms);
%%%%%%%%%%%%%%%%%%%%%%%%
%Vmin= 0*dz;
%Vmax= (Osat-Ohy)*dz;
%Vmax(ms) = Vmax(ms) + (Zbed-ZWT(1))*(Osat-Ohy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%V(1,ms)= 0 + (Osat-Ohy)*(Zbed-ZWT(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(cur_dir)