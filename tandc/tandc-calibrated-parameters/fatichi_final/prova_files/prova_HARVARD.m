%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WORKING LAUNCH PAD HBM  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_directory = cd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
NN= 168120; %%% time Step
%%%%%%%%% Time step 
dt=3600; %%[s] %%% 
dth=1; %%[h]
%%%%%%%%%%%%
ms=14; %%% Soil Layer 
cc = 1; %% Crown area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% METEO INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
id_location = 'HARVARD';
load('E:\SIM_DISC\ARCHIVE\Universit� deposito II\Fluxnet Project\Fluxnet - Data\Harvard_Forest\Data_HARV_run')
NN=length(Date);
%%%%%%%%%%%%%%%%%%%
%N=Latm; 
%%%%%%%%%%%%%%%%%%%
x1=1; % 1562
x2=NN;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%
Ds= esat-ea; 
Date=Date(x1:x2);
Pr=Pr(x1:x2);
Ta=Ta(x1:x2);
Ws=Ws(x1:x2); ea=ea(x1:x2); Pre=Pre(x1:x2); 
Ds=Ds(x1:x2); 
Tdew=Tdew(x1:x2);esat=esat(x1:x2);
SAD1=SAD1(x1:x2); SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2);
SAB2=SAB2(x1:x2);  N=N(x1:x2);
PARB=PARB(x1:x2); PARD = PARD(x1:x2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
t_bef= 0.5; t_aft= 0.5;
%%%%%%%%%%%%%%%%%%%%
Ds(Ds<0)=0;%% [Pa] Vapor Pressure Deficit
%%%%%%%%%%%%%%%%%%%%%%%%% 330-380
%Ca=(393-355)/length(Date)*[1:length(Date)]+ 355; % [ppm]
load('E:\SIM_DISC\DESKTOP_SF\Eco-Hydrology Patterns\CO2_Data\Ca_Data.mat');
d1 = find(abs(Date_CO2-Date(1))<1/36, 1, 'first');d2 = find(abs(Date_CO2-Date(end))<1/36, 1, 'first');
Ca=Ca(d1:d2); 
clear d1 d2 Date_CO2 
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ws(Ws<=0)=0.01;
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM_IC = strcat(current_directory,'\MOD_PARAM_',id_location);
%PARAM_IC = strcat(current_directory,'/MOD_PARAM_',id_location);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Directory = uigetdir('Window','Insert Directory Noname Package') ;
Directory='D:\TeCgam';
cd(Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAIN_FRAME;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('RUN_1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
cd(current_directory);