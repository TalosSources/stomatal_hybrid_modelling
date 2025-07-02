%%%%%%%%%%%% VEGETATION COMPOSITION
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOIL AND HYDROLOGICAL PARAMETER
%%%%%%%%%%%%%%%%%
cur_dir=cd;
cd(Directory)
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
zatm = 30; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%%%%%% VEG. SPECIES  --- Low Grasses -- High Decidous
%%%% LAND COVER PARTITION
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare = 0.0; Ccrown = [1.0];
%%%%%%%%%% SOIL INPUT
Pcla= 0.09;
Psan= 0.65;
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
%%% Measured %%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%%%%5
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
clear Oss_Lp Owp_Lp
Oice = 0;
%%%%%%%%%%%%%
Zs= [0 10 20 50 100 200 300 400 600 800 1000 1200 1500 1750 2000 ]; %% ms+1
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
ExEM = 0.0;
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS VEGETATION
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile
ZR95_H = [1000]; %% [mm]
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
KnitH=0.25; % 0.35; %%% Canopy Nitrogen Decay
KnitL=0.5;
mSl_H = 0.0;%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = 0.0; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_H=[0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[1000]; %%[Pa]
a1_H=[7];
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
%%%%%%%%%%% Matric Potential
%Pss_H = [800]; Pwp_H = [3500]; %%% [kPa]
%Pss_L = [1200]; Pwp_L = [4500]; %%% [kPa]
%%%%%%%% Root Fraction
Psi_sto_00_H =  -0.8; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H =  -2.5 ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H =  -1.0; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H =  -3.0 ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H = 10 ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = 1200;  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = 15.0 ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = 80000;  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = -4.5; %%[MPa]  Water Potential at 50% loss conductivity
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
Sl_H = [0.026]; % 0.024 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [28]; %[kgC/kgN ] Leaf Nitrogen Concentration
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
r_H =  [0.028];  %% [0.066 -0.011] respiration rate at 10° [gC/gN d ]
gR_H=  [0.25];% [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_H=  [ NaN];
aSE_H=  [ 1]; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
dd_max_H=  [ 1/100];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  [ 36/365]; %%20 [Factor of increasing mortality for cold]
Tcold_H =  [2.5]; %%5 [°C] Cold Leaf Shed
drn_H=  [1/1300]; %% turnover root  [1/d]
dsn_H= [1/1000]; % normal transfer rate sapwood [1/d]
age_cr_H= [120]; %% [day] Critical Leaf Age
Bfac_lo_H= [0.95]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [3.5]; %% Mean Temperature for Leaf onset
Tls_H = [NaN]; %% Mean Temperature for Leaf Shed
PAR_th_H= [NaN];
dmg_H= [30]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.01];
Trr_H = [3.0]; %% Translocation rate [gC /m^2 d]
mjDay_H = [250]; %% Maximum Julian day for leaf onset
LDay_min_H =[13.15]; %% Minimum Day duration for leaf onset
LtR_H = [0.3]; %%% Leaf to Root ratio maximum
Mf_H= [1/80]; %% fruit maturation turnover [1/d]
Wm_H= [0] ; % wood turnover coefficient [1/d]
eps_ac_H = [1]; %% Allocation to reserve parameter [0-1]
LDay_cr_H = [11.7]; %%%  Threshold for senescence day light [h]
Klf_H =[1/28]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74]; %% fraction above-ground sapwood and reserve
fbe_H = [0.26]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1]; %% Reference allocation to Fruit and reproduction
[ParEx_H(1)]=Exudation_Parameter(0);
[Mpar_H(1)]=Vegetation_Management_Parameter;
Mpar_H(1).Date_log = datenum(1938,6,15,12,0,0) ; %%%
Mpar_H(1).fract_log = 0.6; %%%
Mpar_H(1).fract_resprout = 0.25;
Mpar_H(1).fract_left=1; %% Fraction_left of leaves and dead leaves
Mpar_H(1).fract_left_fr=1; %% Fraction left of fruits
Mpar_H(1).fract_left_AB = 1; %% Fraction left of harvested/fired wood aboveground
Mpar_H(1).fract_left_BG = 1; %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L = [NaN]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L= [NaN]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
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
Vmax_H = [75]; %62
Vmax_L = [0]; %
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%% Initial Condtion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initial Condtion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%%%%%%%%%%%%%%%
%%% B1 Leaves - Grass  %%% B2 Sapwood  %%% B3 Fine Root  %%% B4 Carbohydrate Reserve
LAI_L(1,:)=[0.0]; %
B_L(1,:,:)= [0 0 0 0 0 0 0 0]; %% 95 58 120 18
Rrootl_L(1,:)= [ 0] ;
PHE_S_L(1,:)=[0];
dflo_L(1,:)=[0];
AgeL_L(1,:)=[0];
e_rel_L(1,:)=[1];
hc_L(1,:) =[0.0]; %% 0.7
SAI_L(1,:) = [0.0]; %% 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAI_H(1,:)=[0];
%B_H(1,:,:)= [ 0 880 450 620 14 0 15 0];
B_H(1,:,:)= [ 0 929 450 663 8.2 11956 0.26 0];
Rrootl_H(1,:)= [4800] ;
PHE_S_H(1,:)=[3];
dflo_H(1,:)=[0];
AgeL_H(1,:)=[0];
e_rel_H(1,:)=[1];
hc_H(1,:) =[23.0];
SAI_H(1,:) = [0.2];
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [1000];
Preserve_H(1,:)= [100];
Kreserve_H(1,:)= [1000];
FNC_H(1,:)=1;
NupI_H(1,:)= [0 0 0];
Nreserve_L(1,:)= [1000];
Preserve_L(1,:)= [100];
Kreserve_L(1,:)= [1000];
FNC_L(1,:)=1;
NupI_L(1,:)= [0 0 0];
RexmyI(1,:)= [0 0 0];
%%%
%%%%%
if OPT_SoilBiogeochemistry == 1
    %%%%% 80.5 years (mean) The site has been regrowing since before 1900
    %%%%% Overstory trees were uprooted by hurricane in 1938, 40–50%
    %%%%% damaged
    %%%% 44.11 yr  --> 14783
    Wm_H= [1/28489]; %% 78 yr //  26157  ; % wood turnover coefficient [1/d]
    %B_H(1,1,6)= 26157;
    %B_H(1,1,6)= 15930 ; % 2007
    %B_H(1,1,6)= 13300 ; % 1991
    %%%
    %%% Pre spin-up
    %B_H(1,:,:)= [ 0 947 450 676 8 11497 0.26 0];
    %%%% After spin-up
    B_H(1,:,:)= [ 0 949 450 678 8.3 12243 0.26 0];
    
    Nreserve_H(1,:)= [19.9];%19.;%
    Preserve_H(1,:)= [2.06];%1.8
    Kreserve_H(1,:)=  [10.2];%10.37; %
    FNC_H(1,:)=[1 ];
    NupI_H(1,:)= [   0.0309   0.0040  0.0112];
    RexmyI(1,:)= [ 0.0398  0.0821        0];
    Nreserve_L(1,:)= [0]; Preserve_L(1,:)= [0]; Kreserve_L(1,:)= [0]; FNC_L(1,:)=[1];
    NavlI(1,:)=[  0.3197    0.1410    0.1353];
    

    
    %%%%%%  After Spin-up
   P(1,:)=  1000*[ 0.125427648115239   0.139662239657644   0.022516778120195   0.902830363616988   0.300943454538995   0.036611460861082,...
   0.144391395501209   0.031753505761769   2.381328863231788   0.680869419502036   5.543308956360717   0.015734528738066,...
   0.011126674058792   0.000062375711503   0.000043029343524   0.000025989879793   0.000071715572540   0.023759780005275,...
   0.093598757940304   0.035282301145315                   0   0.002402860727715   0.003920215994638   0.006294553764415,...
   0.002049439644318   0.699028727459523   0.004574744316396   0.014413739619714   0.001961752896772                   0,...
   0.000191818971091   0.000019222913810   0.000010525538032   0.000095286072772   0.000681273215134   0.000450855946700,...
   0.000228926134117   0.132716895501674   0.001483388729195   0.002343120809116   0.000294262934516                   0,...
   0.000079182087218   0.150000000000000   0.002683103196900   0.015012010597584   0.000000125221416   0.001976839175377,...
   0.003146538554650   0.000293739727691   0.010269280729257   0.000068147054469   0.000065450931968   0.006202387892308,...
   0.503283978894376];
    
    
    %%%% Pre-spin up
   % P(1,:)=  1000*[  0.125262204250535   0.139500594482526   0.022493767059190   0.862753480549643   0.287584493516548   0.036397433848557,...
   %0.138014620549645   0.029814743579083   2.747817884019959   0.673242302721085   5.629335313272005   0.015753772485323,...
   %0.011058785824584   0.000062766208476   0.000043994825255   0.000026152586865   0.000073324708759   0.023910568713284,...
   %0.097222285070500   0.035373159254380                   0   0.002387194871118   0.003906150414526   0.006017448557533,...
   %0.001968846317589   0.723886907391504   0.004603277524930   0.014971263700258   0.001966804754393                   0,...
   %0.000176966351771   0.000017736058337   0.000010651445955   0.000093719487112   0.000595207122099   0.000430415397024,...
   %0.000198051378413   0.134565736283261   0.001492656126511   0.002433606475366   0.000295020713159                   0,...
   %0.000033018333743   0.150000000000000   0.000748913245098   0.015006365590690   0.000000123017705   0.001988675277407,...
   %0.003009264056597   0.000285677518164   0.010360348353737   0.000075275850610   0.000072288385688   0.007951685353006,...
   %0.503263370750018];
    
    
    Upl=0.01; %% Soil production [mm/yr]
    HIST=0;
    %[B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl,HIST,FertN,FertP,FertK,DepN,DepP,DepK);
    [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl,HIST);
    B_IO.SC_par=[1 1 1 1];
    PHs=4.3;  %
end
%%%
TBio_L=0;  %%[ton DM / ha ]
TBio_H=200; %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=10;  %% [mm/ m2 PFT];
Vl_H=10;  %% [mm/ m2 PFT];
Vx_L=0;   %% [mm/ m2 PFT];
Vl_L=0;   %% [mm/ m2 PFT];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+2;
Tdamp(1)=Ta(1);
Tdp(1,:)= 4.5*ones(1,ms);
TdpI_H(1,:)=4;
TdpI_L(1,:)=4;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
cd(cur_dir)
%%%%%%%%%%%%%%%%%