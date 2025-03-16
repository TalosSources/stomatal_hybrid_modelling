%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WORKING LAUNCH PAD HBM  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INPUT MANAGER 
current_directory = cd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
NN=8760;% 276096;% %%% time Step
%%%%%%%%% Time step 
dt=3600; %%[s] %%% 
dth=1; %%[h]
%%%%%%%%%%%%
ms=16; %%% Soil Layer 
cc = 1; %% Crown area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% METEO INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_location = 'ZURICH_SMA';
load('C:\Users\tlian\Desktop\Soil\T&C\T&C_tutorial\T&C_file_model\Inputs\Data_Run_Zurich_Fluntern.mat')
Date=D; clear D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
x1=1;
x2=8760;% 276096;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%5
%% Meteorological inputs
Date=Date(x1:x2);
Pr=Pr(x1:x2); % Precipitation [mm/h]
Ta=Ta(x1:x2); % Air Temperature 
Ws=Ws(x1:x2); ea=ea(x1:x2);  SAD1=SAD1(x1:x2); % wind speed [m/s] Vapor Pressure  first band Diffuse radiation
SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2); Pre=Pre(x1:x2); %Atmospheric Pressure
SAB2=SAB2(x1:x2); N=N(x1:x2); Tdew=Tdew(x1:x2);esat=esat(x1:x2); % N:Cloud Cover or 
% Longwave Incoming Radiation  Dew Point Temperature 
PARB=PARB(x1:x2); PARD = PARD(x1:x2); % PAR(Photosynthetically Active Radiation) radiation Direct  PAR radiation Diffuse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
t_bef= -0.67; t_aft= 1.67; %Integration interval for solar variables – Hours or fraction, based on monitor data
%%%%%%%%%%%%%%%%%%%%
Ds=esat-ea; %% [Pa] Vapor Pressure Deficit
Ds(Ds<0)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%
load('C:\Users\tlian\Desktop\Soil\T&C\T&C_tutorial\T&C_file_model\Inputs\Ca_Data.mat');
d1 = find(abs(Date_CO2-Date(1))<1/36);d2 = find(abs(Date_CO2-Date(end))<1/36);
Ca=Ca(d1:d2); %CO2 atmospheric concentration
clear d1 d2 Date_CO2 
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ws(Ws<=0)=0.01;
%%dt,Pr(i),Ta(i),Ws(i),ea(i),Pre(i),Rdir(i),Rdif(i),N(i),z,Tdew(i),esat(i),.
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM_IC = strcat(current_directory,'\MOD_PARAM_',id_location);
%PARAM_IC = strcat(current_directory,'/PARAMETER_PLOT_FILES/MOD_PARAM_',id_location);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Directory = uigetdir('Window','Insert Directory Noname Package') ;
Directory='C:\Users\tlian\Desktop\Soil\T&C\T&C_tutorial\T&C_file_model\T&C_CODE';
cd(Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAIN_FRAME ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(current_directory);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_to_save = strcat('C:\Users\tlian\Desktop\Soil\T&C\T&C_tutorial\T&C_file_model\results\Ris_',id_location,'.mat');
save(file_to_save);
%%%%%%%%%%%%%%%%%%%%%%%%