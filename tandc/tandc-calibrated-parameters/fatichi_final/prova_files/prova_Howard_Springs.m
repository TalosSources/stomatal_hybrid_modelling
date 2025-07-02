%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WORKING LAUNCH PAD HBM  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INPUT MANAGER 
current_directory = cd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
NN= 122712;%%% time Step
%%%%%%%%% Time step 
dt=3600; %%[s] %%% 
dth=1; %%[h]
%%%%%%%%%%%%
ms=27; %%% Soil Layer 
cc = 2; %% Crown area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% METEO INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_location = 'Howard_Springs';
load('E:\SIM_DISC\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\New_Sites\Data_Howard_Springs_Run.mat')
NN=length(Date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('E:\SIM_DISC\desk\Savannah\Data_Howard_Springs_detrended.mat')
%Pr=Pr_noIV; 
%load('E:\SIM_DISC\desk\Savannah\Data_Howard_Springs_novariability.mat')
%load('E:\SIM_DISC\desk\Savannah\Data_Howard_Springs_novariability_Pr.mat')
% Ca= Ca_rep ;
% SAB1 = SAB1_rep ;
% SAB2=SAB2_rep;
% SAD1=SAD1_rep;
% SAD2=SAD2_rep;
% PARB=PARB_rep;
% PARD=PARD_rep;
% Ws=Ws_rep;
% U=U_rep;
% ea = ea_rep;
% esat = esat_rep; 
% Tdew = Tdew_rep; 
%Pr=Pr_rep;
%Ta=Ta_rep;
%Pr=Pr_rep_exc_270_300; 
%%%%%%%%%%%%%%%%%%%
x1=1;
x2=NN;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%5
Date=Date(x1:x2);
Pr=Pr(x1:x2);
Ta=Ta(x1:x2);
Ws=Ws(x1:x2); ea=ea(x1:x2);  SAD1=SAD1(x1:x2);
SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2);
SAB2=SAB2(x1:x2); N=N(x1:x2); Tdew=Tdew(x1:x2);esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD = PARD(x1:x2); Pre=Pre(x1:x2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
t_bef=0.50; t_aft=  0.50;  
%%%%%%%%%%%%%%%%%%%% 
%Pre=ones(1,NN)*1000; %% [mbar]
Ds=esat-ea; %% [Pa] Vapor Pressure Deficit
Ds(Ds<0)=0;
%%%%%%%%%%%%%%%%%%%%%%%%% 330-380
load('E:\SIM_DISC\DESKTOP_SF\Eco-Hydrology Patterns\CO2_Data\Ca_Data.mat');
d1 = find(abs(Date_CO2-Date(1))<1/36);d2 = find(abs(Date_CO2-Date(end))<1/36);
Ca=Ca(d1:d2); 
clear d1 d2 Date_CO2 
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ws(Ws<=0)=0.01;
%%dt,Pr(i),Ta(i),Ws(i),ea(i),Pre(i),Rdir(i),Rdif(i),N(i),z,Tdew(i),esat(i),.
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
MAIN_FRAME ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('RUN_1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
cd(current_directory);