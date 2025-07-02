%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOIL AND HYDROLOGICAL PARAMETER
cur_dir=cd;
cd(Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_dis=_a_dis_ ;
pow_dis=_pow_dis_;
%%%%%%%%%%%%%%%%%%%%
fpr=_fpr_;
SvF=_SvF_; %% Sky View Factor
SN=_SN_; %% Stream Identifier
Slo_top=_Slo_top_;  %% [fraction dy/dx]
Slo_pot=zeros(1,ms); %% [fraction dy/dx]
Asur=1./cos(atan(Slo_top)); %% Real Area/Projected Area [m^2/m^2]
Ared=_Ared_;
aR=_aR_; %%% anisotropy ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsize=_cellsize_; %%[m^2];
aTop=1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
Kbot=_Kbot_; %% [mm/h] Conductivity at the bedrock layer
Krock=_Krock_; %% [mm/h] Conductivity of Fractured Rock
zatm=_zatm_; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%%% LAND COVER PARTITION
Cwat=_Cwat_; 
Curb=_Curb_; 
Crock=_Crock_;
Cbare=_Cbare_; 
Ccrown=_Ccrown_;
%%%%%%%%%% SOIL INPUT
Pcla=_Pcla_;
Psan=_Psan_;
Porg=_Porg_; %
Color_Class=_Color_Class_;
%%%%%%%%%%%%%%%%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%%%%%%%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s=lan_s*ones(1,ms);
cv_s=cv_s*ones(1,ms);
%%%
%%%%
SPAR=_SPAR_; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%
%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe=Pe*ones(1,ms);
O33=O33*ones(1,ms);
alpVG=alpVG*ones(1,ms); %% [1/mm]
nVG=nVG*ones(1,ms); %% [-]
Ks_Zs=Ks*ones(1,ms); %%[mm/h]
%%%%%%%%%%% Matric Potential
Kfc=_Kfc_; %% [mm/h]
Phy=_Phy_; %% [kPa]
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
clear Oss_Lp Owp_Lp
Oice = _Oice_;
%%%%%%%%%%%%%%%%%%%
Zs=_Zs_; %%% [ms+1]
Zdes=_Zdes_;
Zinf=_Zinf_;
Zbio = _Zbio_;
if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end
[EvL_Zs]=Evaporation_layers(Zs,Zdes); %%% Evaporation Layer fraction
[Inf_Zs]=Evaporation_layers(Zs,Zinf); %%% Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio); %%% Infiltration Depth Layer fraction
dz=diff(Zs); %%%% [mm]  Thickness of the Layers
Dz=zeros(1,ms);
for i = 1:ms
    if i>1
        Dz(i)= (dz(i)+ dz(i-1))/2; %%% Delta Depth Between Middle Layer  [mm]
    else
        Dz(i)=dz(1)/2; %%% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end
%%%%%%%%%%%%%%%%% OTHER PARAMETER
In_max_urb=_In_max_urb_; In_max_rock=_In_max_rock_; %% [mm]
%%%%%%%%%%%%% SNOW PARAMETER
TminS=_TminS_;%% Threshold temperature snow
TmaxS=_TmaxS_;%% Threshold temperature snow
ros_max1=_ros_max1_; %520 600; %%% [kg/m^3]
ros_max2=_ros_max2_; %320 450; %%% [kg/m^3]
Th_Pr_sno=_Th_Pr_sno_; %%% [mm/day] Threshold Intensity of snow to consider a New SnowFall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ICE Parameter
Ice_wc_sp=_Ice_wc_sp_; %% [-] Specific Maximum water content ice
ros_Ice_thr=_ros_Ice_thr_ ; %% [kg/m^3] Density Thrshold to transform snow into ice
Aice=_Aice_; %% [-] Ice albedo
WatFreez_Th=_WatFreez_Th_; %% [°C] Threshold for freezing lake water
dz_ice=_dz_ice_; %% [mm / h] Water Freezing Layer progression without snow-layer
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
ExEM=_ExEM_;
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS VEGETATION
%%% cc -- number of crown area
CASE_ROOT=_CASE_ROOT_;  %%% Type of Root Profile
ZR95_H=_ZR95_H_; %% [mm]
ZR95_L=_ZR95_L_; %% [mm]
ZR50_H=_ZR50_H_;
ZR50_L=_ZR50_L_;
ZRmax_H=_ZRmax_H_;
ZRmax_L=_ZRmax_L_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kct=_Kct_; %%% Factor Vegetation Cover --- for throughfall
gcI=_gcI_; %%% [1/mm]
KcI=_KcI_; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In=_Sp_SN_In_; %% [mm/LAI]
Sp_LAI_H_In=_Sp_LAI_H_In_; %%[mm/LAI]
Sp_LAI_L_In=_Sp_LAI_L_In_; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H=_d_leaf_H_; %%[cm]
d_leaf_L=_d_leaf_L_;  %% [cm]
%%%%%%%% Biochemical parameter
KnitH=_KnitH_; %%% Canopy Nitrogen Decay
KnitL=_KnitL_;
mSl_H=_mSl_H_;%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L=_mSl_L_; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Photosynthesis Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_H=_FI_H_;% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=_Do_H_; %%[Pa]
a1_H=_a1_H_; % influences transpiration - latent heat
go_H=_go_H_;% [mol / s m^2] minimum Stomatal Conductance
CT_H=_CT_H_; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H=_DSE_H_;  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H=_Ha_H_; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=_gmes_H_; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=_rjv_H_; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FI_L=_FI_L_;% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=_Do_L_; %%[Pa]
a1_L=_a1_L_;
go_L=_go_L_;% % [mol / s m^2] minimum Stomatal Conductance
CT_L=_CT_L_;  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L=_DSE_L_;  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L=_Ha_L_; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=_gmes_L_; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L=_rjv_L_; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pss_H = [1000]; Pwp_H = [2500]; %%% [kPa]
%Pss_L=50; Pwp_L=300; %%% [kPa]
Psi_sto_00_H=_Psi_sto_00_H_; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H=_Psi_sto_50_H_ ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H=_PsiL00_H_; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H=_PsiL50_H_;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H=_Kleaf_max_H_ ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H=_Cl_H_;  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H=_Axyl_H_ ; %% [cm^2 stem /m^2 PFT]
Kx_max_H=_Kx_max_H_ ;  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H=_PsiX50_H_; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H=_Cx_H_; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata
Psi_sto_00_L=_Psi_sto_00_L_;%  %% [MPa]  Water Potential at PLCs loss conductivity
Psi_sto_50_L=_Psi_sto_50_L_;%  %% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L=_PsiL00_L_; %% [MPa]  Water Potential at PLCs% loss conductivity
PsiL50_L=_PsiL50_L_;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L=_Kleaf_max_L_ ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L=_Cl_L_;  %%%  [500 - 3000]%  [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L=_Axyl_L_ ; %% [cm^2 stem /m^2 PFT]
Kx_max_L=_Kx_max_L_;  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L=_PsiX50_L_; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L=_Cx_L_; %%% [kg / m^3 sapwood MPa]
%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);

%%%% Growth Parameters
PsiG50_H=_PsiG50_H_;  %%[MPa]
PsiG99_H=_PsiG99_H_;  %%[MPa]
gcoef_H=_gcoef_H_; % [gC/m2 day]
%%------
PsiG50_L=_PsiG50_L_;
PsiG99_L=_PsiG99_L_;
gcoef_L=_gcoef_L_; % [gC/m2 day]
%%%%%%%% Vegetation Optical Parameter
_Comment_PFT_opt_H_1_[PFT_opt_H(1)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_H_1_);
_Comment_PFT_opt_H_2_[PFT_opt_H(2)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_H_2_);
_Comment_PFT_opt_H_3_[PFT_opt_H(3)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_H_3_);
_Comment_PFT_opt_H_4_[PFT_opt_H(4)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_H_4_);
_Comment_PFT_opt_H_5_[PFT_opt_H(5)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_H_5_);
_Comment_PFT_opt_H_6_[PFT_opt_H(6)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_H_6_);

_Comment_PFT_opt_L_1_[PFT_opt_L(1)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_L_1_);
_Comment_PFT_opt_L_2_[PFT_opt_L(2)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_L_2_);
_Comment_PFT_opt_L_3_[PFT_opt_L(3)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_L_3_);
_Comment_PFT_opt_L_4_[PFT_opt_L(4)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_L_4_);
_Comment_PFT_opt_L_5_[PFT_opt_L(5)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_L_5_);
_Comment_PFT_opt_L_6_[PFT_opt_L(6)]=Veg_Optical_Parameter(_Veg_Optical_Parameter_L_6_);

OM_H=_OM_H_;
OM_L=_OM_L_;

Sllit=_Sllit_ ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H=_Sl_H_ ; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H=_Nl_H_; %[gC/gN ] Leaf Carbon-Nitrogen ratio
_Comment_Stoich_H_1_[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
_Comment_Stoich_H_2_[Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
_Comment_Stoich_H_3_[Stoich_H(3)]=Veg_Stoichiometric_Parameter(Nl_H(3));
_Comment_Stoich_H_4_[Stoich_H(4)]=Veg_Stoichiometric_Parameter(Nl_H(4));
_Comment_Stoich_H_5_[Stoich_H(5)]=Veg_Stoichiometric_Parameter(Nl_H(5));
_Comment_Stoich_H_6_[Stoich_H(6)]=Veg_Stoichiometric_Parameter(Nl_H(6));
%PLNR_H = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_H=_r_H_;  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_H=_gR_H_; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_H= [2.0 10.0];
aSE_H=_aSE_H_ ; % 0; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Trop Eve. For. 
dd_max_H=_dd_max_H_; %%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H=_dc_C_H_; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H=_Tcold_H_; %% [°C] Cold Leaf Shed
drn_H=_drn_H_; %% turnover root  [1/d]
dsn_H=_dsn_H_; % normal transfer rate sapwood [1/d]
age_cr_H=_age_cr_H_;% 365; %% [day] Critical Leaf Age
Bfac_lo_H=_Bfac_lo_H_; %% Leaf Onset Water Stress
Bfac_ls_H=_Bfac_ls_H_; %% Leaf Shed Water Stress [0-1]
Tlo_H=_Tlo_H_; %% Mean Temperature for Leaf onset
Tls_H=_Tlo_H_; %% Mean Temperature for Leaf Shed
PAR_th_H=_PAR_th_H_; 
dmg_H=_dmg_H_ ;% 30; %%%  Day of Max Growth
LAI_min_H=_LAI_min_H_;
Trr_H=_Trr_H_; % 1; %% Translocation rate [gC /m^2 d]
mjDay_H=_mjDay_H_; %% 180 %% Maximum Julian day for leaf onset
LDay_min_H=_LDay_min_H_; %% Minimum Day duration for leaf onset
LtR_H=_LtR_H_; %%% Leaf to Root ratio maximum
Mf_H=_Mf_H_; %% fruit maturation turnover [1/d]
Wm_H=_Wm_H_; % wood turnover coefficient [1/d]
eps_ac_H=_eps_ac_H_; %% Allocation to reserve parameter [0-1]
LDay_cr_H=_LDay_cr_H_; %%%  Threshold for senescence day light [h]
Klf_H=_Klf_H_; %% Dead Leaves fall turnover [1/d]
fab_H=_fab_H_; %% fraction above-ground sapwood and reserve
fbe_H=_fbe_H_; %% fraction below-ground sapwood and reserve
ff_r_H=_fr_H_; %% Reference allocation to Fruit and reproduction
_Comment_ParEx_H_1_[ParEx_H(1)]=Exudation_Parameter(_Exudation_Parameter_H_1_);  
_Comment_ParEx_H_2_[ParEx_H(2)]=Exudation_Parameter(_Exudation_Parameter_H_2_);  
_Comment_ParEx_H_3_[ParEx_H(3)]=Exudation_Parameter(_Exudation_Parameter_H_3_);  
_Comment_ParEx_H_4_[ParEx_H(4)]=Exudation_Parameter(_Exudation_Parameter_H_4_);  
_Comment_ParEx_H_5_[ParEx_H(5)]=Exudation_Parameter(_Exudation_Parameter_H_5_);  
_Comment_ParEx_H_6_[ParEx_H(6)]=Exudation_Parameter(_Exudation_Parameter_H_6_);  
_Comment_Mpar_H_1_[Mpar_H(1)]=Vegetation_Management_Parameter;
_Comment_Mpar_H_2_[Mpar_H(2)]=Vegetation_Management_Parameter;
_Comment_Mpar_H_3_[Mpar_H(3)]=Vegetation_Management_Parameter;
_Comment_Mpar_H_4_[Mpar_H(4)]=Vegetation_Management_Parameter;
_Comment_Mpar_H_5_[Mpar_H(5)]=Vegetation_Management_Parameter;
_Comment_Mpar_H_6_[Mpar_H(6)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L=_Sl_L_; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L=_Nl_L_; %[kgC/kgN ] Leaf Nitrogen Concentration
_Comment_Stoich_L_1_[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
_Comment_Stoich_L_2_[Stoich_L(2)]=Veg_Stoichiometric_Parameter(Nl_L(2));
_Comment_Stoich_L_3_[Stoich_L(3)]=Veg_Stoichiometric_Parameter(Nl_L(3));
_Comment_Stoich_L_4_[Stoich_L(4)]=Veg_Stoichiometric_Parameter(Nl_L(4));
_Comment_Stoich_L_5_[Stoich_L(5)]=Veg_Stoichiometric_Parameter(Nl_L(5));
_Comment_Stoich_L_6_[Stoich_L(6)]=Veg_Stoichiometric_Parameter(Nl_L(6));
%PLNR_L = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_L=_r_L_;  %% [0.066 -0.011] respiration rate at 10° [gC/gN d ]
gR_L=_gR_L_;% [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_L=  [ NaN NaN];
aSE_L=_aSE_L_; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
dd_max_L=_dd_max_L_;%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L=_dc_C_L_; %% [Factor of increasing mortality for cold]
Tcold_L=_Tcold_L_; %% [°C] Cold Leaf Shed
drn_L=_drn_L_; %% senescence rate leaf  [1/d]
dsn_L=_dsn_L_; % [ normal transfer rate sapwood [1/d] ]
age_cr_L=_age_cr_L_; %% [day] Critical Leaf Age
Bfac_lo_L=_Bfac_lo_L_; %% PHENOLOGY Leaf Onset Water Stress [0-1]
Bfac_ls_L=_Bfac_ls_L_; %% PHENOLOGY Leaf Shed Water Stress [0-1]
Tlo_L=_Tlo_L_; %% Mean Temperature for Leaf onset
Tls_L=_Tls_L_; %% Mean Temperature for Leaf Shed
PAR_th_L=_Tls_L_; 
dmg_L=_dmg_L_; %%% Tree 30 Grasses Day of Max Growth
LAI_min_L=_LAI_min_L_;
Trr_L=_Trr_L_; %% Translocation rate [gC /m^2 d]
mjDay_L=_mjDay_L_; %% Maximum Julian day for leaf onset
LDay_min_L=_LDay_min_L_; %% Minimum Day duration for leaf onset
LtR_L=_LtR_L_; %%% Leaf to Root ratio maximum
Mf_L=_Mf_L_; %% fruit maturation turnover [1/d]
Wm_L=_Wm_L_; % wood turnover coefficient [1/d]
eps_ac_L=_eps_ac_L_; %% Allocation to reserve parameter [0-1]
LDay_cr_L=_LDay_cr_L_ ; % 12.75; %%%  Threshold for senescence day light [h]
Klf_L=_Klf_L_; %% Dead Leaves fall turnover [1/d]
fab_L=_fab_L_; %% fraction above-ground sapwood and reserve
fbe_L=_fbe_L_; %% fraction below-ground sapwood and reserve
ff_r_L=_fr_L_;
_Comment_ParEx_L_1_[ParEx_L(1)]=Exudation_Parameter(_Exudation_Parameter_L_1_);  
_Comment_ParEx_L_2_[ParEx_L(2)]=Exudation_Parameter(_Exudation_Parameter_L_2_);  
_Comment_ParEx_L_3_[ParEx_L(3)]=Exudation_Parameter(_Exudation_Parameter_L_3_);  
_Comment_ParEx_L_4_[ParEx_L(4)]=Exudation_Parameter(_Exudation_Parameter_L_4_);  
_Comment_ParEx_L_5_[ParEx_L(5)]=Exudation_Parameter(_Exudation_Parameter_L_5_);  
_Comment_ParEx_L_6_[ParEx_L(6)]=Exudation_Parameter(_Exudation_Parameter_L_6_);  
_Comment_Mpar_L_1_[Mpar_L(1)]=Vegetation_Management_Parameter;
_Comment_Mpar_L_2_[Mpar_L(2)]=Vegetation_Management_Parameter;
_Comment_Mpar_L_3_[Mpar_L(3)]=Vegetation_Management_Parameter;
_Comment_Mpar_L_4_[Mpar_L(4)]=Vegetation_Management_Parameter;
_Comment_Mpar_L_5_[Mpar_L(5)]=Vegetation_Management_Parameter;
_Comment_Mpar_L_6_[Mpar_L(6)]=Vegetation_Management_Parameter;
%Mpar_L(1).jDay_cut =[125  156  186  217  247  278];
%Mpar_L(1).LAI_cut =[1.68]; %% LAI of grass after cut %% 2.65 Height of 15cm more or less
%%%%%%%%%%%%% PRODUCTIVITY
%%% Maximum Rubisco Capacity
%[Vmax_H]=Maximum_Rubisco_Capacity(Sl_H,PLNR_H,Nl_H);
%[Vmax_L]=Maximum_Rubisco_Capacity(Sl_L,PLNR_L,Nl_L);
Vmax_H=_Vmax_H_; %55
Vmax_L=_Vmax_L_; %35
%[Amax_H]= MAX_PHOTOSYNTESIS(Vmax_H,Ca,CT_H,Tup_H,Tlow_H,FI_H,Oa); %% 17
%[Amax_L]= MAX_PHOTOSYNTESIS(Vmax_L,Ca,CT_L,Tup_L,Tlow_L,FI_L,Oa); %% 13
%%%%%%%%%%%%%%%%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')

%%%% Initial Condtion
%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:)=Ca(1); % %% [umolCO2/mol]
Ci_sunH(1,:)=Ca(1); %% [umolCO2/mol]
Ci_shdL(1,:)=Ca(1); % %% [umolCO2/mol]
Ci_shdH(1,:)=Ca(1); %% [umolCO2/mol]
%%%%%%%%%%%%%%%%
LAI_H(1,:)=_Init_LAI_H_; %
B_H(1,:,:)=_Init_B_H_; %%
Rrootl_H(1,:)=_Init_Rrootl_H_;
PHE_S_H(1,:)=_Init_PHE_S_H_;
dflo_H(1,:)=_Init_dflo_H_;
AgeL_H(1,:)=_Init_AgeL_H_;
e_rel_H(1,:)=_Init_e_rel_H_;
hc_H(1,:)=_Init_hc_H_; %% 0.7
SAI_H(1,:)=_Init_SAI_H_; %% 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LAI_L(1,:)=_Init_LAI_L_;
B_L(1,:,:)=_Init_B_L_;
Rrootl_L(1,:)=_Init_Rrootl_L_;
PHE_S_L(1,:)=_Init_PHE_S_L_;
dflo_L(1,:)=_Init_dflo_L_;
AgeL_L(1,:)=_Init_AgeL_L_;
e_rel_L(1,:)=_Init_e_rel_L_;
hc_L(1,:)=_Init_hc_L_;
SAI_L(1,:)=_Init_SAI_L_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=_Init_BLit_;  %% [kg DM /m2 PFT] Litter Biomass
%%%%
PARI_H(1,:,:)=_Init_PARI_H_;
NBLI_H(1,:)=_Init_NBLI_H_;
PARI_L(1,:,:)=_Init_PARI_L_; 
NBLI_L(1,:)=_Init_NBLI_L_;
%%%
Nreserve_H(1,:)=_Init_Nreserve_H_;
Preserve_H(1,:)=_Init_Preserve_H_;
Kreserve_H(1,:)=_Init_Kreserve_H_;
FNC_H(1,:)=_Init_FNC_H_;
NupI_H(1,:)=_Init_NupI_H_;
Nreserve_L(1,:)=_Init_Nreserve_L_;
Preserve_L(1,:)=_Init_Preserve_L_;
Kreserve_L(1,:)=_Init_Kreserve_L_;
FNC_L(1,:)=_Init_FNC_L_;
NupI_L(1,:)=_Init_NupI_L_;
RexmyI(1,:)=_Init_RexmyI_;
%%%%
%%%
if OPT_SoilBiogeochemistry == 1
    %%%%%
    Nreserve_L(1,:)= 7.5; %
    Preserve_L(1,:)= 0.55 ; %
    Kreserve_L(1,:)= 2.8 ; %
    FNC_L(1,:)=1;
    NupI_L(1,:)= [ 0.10866   0.007844   0.041600];
    RexmyI(1,:)= [ 0.04685   0.12003              0];
    NavlI(1,:)=[   0.491   0.0363   0.1653];
    %%%%
    %%%
    P(1,:)=      1000*[ 0.125058323623613   0.331492841956637   0.055153164737343                   0                   0   0.106545776397086,...
   0.130152646652338   0.014461405183593   0.679083359004130   0.966530132450958   6.682770254265284   0.008562021215738,...
   0.006309642813074   0.000068954335702   0.000043373075760   0.000028730973209   0.000072288459601   0.043379971628909,...
   0.142594653474930   0.051308356465522                   0   0.002269771399924   0.010807463358596                   0,...
   0.004809032005655   0.946025338170192   0.008369508208843   0.021979657773539   0.002852827441936                   0,...
   0.000309843383482   0.000031633522226   0.000026495680780   0.000226977139993   0.000780834292220                   0,...
   0.000346262009287   0.137514873573224   0.002709812306318   0.003570543222058   0.000427924116290                   0,...
   0.000022069475017   0.150000000000000   0.000834499816223   0.015010533133807   0.000000222487831   0.005435212585202,...
                   0   0.000558301361478   0.020245150099706   0.000040438105769   0.000033986846738   0.004410679675017,...
   0.503240898696182]; 

    %%% SON SOP  946 136.7
    %%%
    %DepN= 0.00136; %
    %DepP=1.578e-005; %  gP/m2 day
    %DepK=1.27e-004;
    %FertN=0*ones(1,366);
    %FertP=0*ones(1,366);
    %FertK=0*ones(1,366);
    Upl=0.01; %% Soil production [mm/yr]
    HIST=0;
    %[B_IO]=Biogeochemistry_IO(Zs(ms+1)*1000,Lat,Lon,Upl,HIST,FertN,FertP,FertK,DepN,DepP,DepK);
    [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl);
    B_IO.SC_par=[1 1 1 1];
    PHs=5; 
end
%%%
TBio_L=_TBio_L_;  %%[ton DM / ha ]
TBio_H=_TBio_H_;  %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=_Vx_H_;  %% [mm/ m2 PFT];
Vl_H=_Vl_H_;  %% [mm/ m2 PFT];
Vx_L=_Vx_L_;   %% [mm/ m2 PFT];
Vl_L=_Vl_L_;   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%%%%%%% Initial Conditions
SWE(1)=_Init_SWE_; %% [mm]
SND(1)=_Init_SND_;
Ts(1)=Ta(1)+2;
Tdamp(1)=Ta(1);
Tdp(1,:)=Ta(1)*ones(1,ms);
Sdp(1,:)=(0.001*dz).*cv_s.*(Tdp(1,:)); %% [J °C/m^2 K]
TdpI_H(1,:)=_TdpI_H_;
TdpI_L(1,:)=_TdpI_L_;
%%% Snow_alb = soil_alb initial
snow_alb.dir_vis=_Init_snow_alb_dir_vis_;
snow_alb.dif_vis=_Init_snow_alb_dif_vis_;
snow_alb.dir_nir=_Init_snow_alb_dir_nir_;
snow_alb.dif_nir=_Init_snow_alb_dif_nir_;
In_L(1,:)=_In_L_; In_H(1,:)=_In_H_;
In_urb(1)=_In_urb_; In_rock(1)=_In_rock_;
In_Litter(1)=_In_Litter_;
SP_wc(1)=_Init_SP_wc_; %%[mm]
In_SWE(1)=_In_SWE_;
ros(1)=_Init_ros_;
t_sls(1)=_Init_t_sls_;
e_sno(1)=_Init_e_sno_;
tau_sno(1)=_Init_tau_sno_;
EK(1)=_Init_EK_;
WAT(1)=_Init_WAT_;
ICE(1)=_Init_ICE_;
IP_wc(1)=_Init_IP_wc_;
ICE_D(1)=_Init_ICE_D_;
FROCK(1)=_Init_FROCK_;
Ws_under(1)=_Init_Ws_under_;
%%%%%%%%%%%%%% Volume [mm]
%%%%%%%%%%%%%%%%%%%%%
%ZWT(1)= 1999; %% [mm]
O(1,:)=Ofc;
%dz(ms)=ZWT(1)-Zs(ms);
%%%%%%%%%%%%%%%%%%%%%%%%
%Vmin= 0*dz;
%Vmax= (Osat-Ohy)*dz;
%Vmax(ms) = Vmax(ms) + (Zbed-ZWT(1))*(Osat-Ohy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V(1,:)=(O(1,:)-Ohy).*dz;
%V(1,ms)= 0 + (Osat-Ohy)*(Zbed-ZWT(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(cur_dir)
