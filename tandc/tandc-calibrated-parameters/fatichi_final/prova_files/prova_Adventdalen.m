%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WORKING LAUNCH PAD HBM  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_directory = cd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
NN= 35064;  %%% time Step 
%%%
dt=3600; %%[s] %%% 
dth=1; %%[h]
%%%%%%%%%%%%
ms=21; %%% Soil Layer 
cc = 1; %% Crown area 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% METEO INPUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_location = 'Adventdalen'; 
load('C:\DD\ARCHIVE\Universit� deposito II\Fluxnet Project\Fluxnet - Data\FLX_SJ-Adv_FLUXNET2015_SUBSET_2011-2014_1-4\Data_Adventdalen_run.mat')
NN = length(Date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%Wlev = NaN*Pr;  %%% [mm] Wetland Water level series // determined by runoff 
N=Latm; 
%%%%%%%%%%%%%%%%%%%
x1=1; 
x2=NN;  
%%%%%%%%%%%%%%%%%%%%%%%%
Date=Date(x1:x2); 
Pr=Pr(x1:x2);
Ta=Ta(x1:x2);
Ws=Ws(x1:x2); ea=ea(x1:x2); Pre=Pre(x1:x2); SAD1=SAD1(x1:x2); 
SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2); 
SAB2=SAB2(x1:x2); N=N(x1:x2);Tdew=Tdew(x1:x2);esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD = PARD(x1:x2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
t_bef= 0.0; t_aft= 1.0;
%%%%%%%%%%%%%%%%%%%%
Ds=esat-ea; %% [Pa] Vapor Pressure Deficit 
Ds(Ds<0)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('C:\DD\DESKTOP_SF\Eco-Hydrology Patterns\CO2_Data\Ca_Data.mat');
d1 = find(abs(Date_CO2-Date(1))<1/36, 1, 'first');d2 = find(abs(Date_CO2-Date(end))<1/36, 1, 'first');
Ca=Ca(d1:d2); 
clear d1 d2 Date_CO2 
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
Ws(Ws<=0)=0.01; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%dt,Pr(i),Ta(i),Ws(i),ea(i),Pre(i),Rdir(i),Rdif(i),N(i),z,Tdew(i),esat(i),.
[YE,MO,DA,HO,MI,SE] = datevec(Date); 
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO; 
clear YE MO DA HO MI SE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
PARAM_IC = strcat(current_directory,'\MOD_PARAM_',id_location); 
%PARAM_IC = strcat(current_directory,'/MOD_PARAM_',id_location); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Directory = uigetdir('Window','Insert Directory Noname Package') ;
Directory='C:\DD\DESKTOP_SF\BILANCIO IDROLOGICO\TeCgam\TeC_Code';
cd(Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAIN_FRAME ; 
save('RUN_1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(current_directory); 