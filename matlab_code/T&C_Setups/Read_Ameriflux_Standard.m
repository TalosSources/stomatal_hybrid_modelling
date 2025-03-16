%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Run Section by Section!!!
%--------------------------------
%The following code is reading Ameriflux half hourly (HH) data %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%Site description
Lon          = -110.8661;   
Lat          = 31.8214;
DeltaGMT     = -7;
Zbas         = 1120;
Location_Name ='Santa_Rita_Mesquite'; 
%%%%%%%%%%%%%%%%%


cur_dir=cd; 

file_dir = ('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\AMF_US-SRM_BASE-BADM_16-5');
file_name = 'AMF_US-SRM_BASE_HH_16-5.csv';
file_name_less = file_name(1:end-4);


DATA=csvread(strcat(file_dir,'\',file_name),3,0);
DATA(DATA==-9999)=NaN;

% [p,c] =xlsread(strcat(file_dir,'\',file_name),file_name_less,'A3:A3');
% c=char(c);
% 
% Var_List = strsplit(c,',');
% clear pc c

[p,c] =xlsread(strcat(file_dir,'\',file_name),file_name_less,'A3:AL3');
%c=char(c);

Var_List=c; 
%Var_List = strsplit(c,',');
clear pc c

%TIMESTAMP_START,TIMESTAMP_END,USTAR,TA,WD,WS,NEE_PI,FC,H,LE,TS_PI_1,TS_PI_2,RH,PA,CO2_PI_1,CO2_PI_2,VPD_PI,SWC_PI_1,SWC_PI_2,NETRAD,PPFD_IN,SW_IN,SW_OUT,LW_IN,LW_OUT,
%H2O,RECO_PI,PPFD_DIF,CO2,SC,NEE_PI_F,
%RECO_PI_F,GPP_PI,GPP_PI_F,TS_1_1_1,TS_1_2_1,TS_1_3_1,TS_1_4_1,TS_1_5_1,TS_1_6_1,TS_1_7_1,TS_1_8_1,SWC_1_1_1,SWC_1_2_1,SWC_1_3_1,SWC_1_4_1,SWC_1_5_1,SWC_1_6_1,SWC_1_7_1,SWC_1_8_1

%TIMESTAMP_START,TIMESTAMP_END,NETRAD,G,LE,H,FC,USTAR,CO2,H2O,T_SONIC,PA,TA_1_2_1,
%RH_1_2_1,WS_1_2_1,WD_1_2_1,TA_1_1_1,RH_1_1_1,WS_1_1_1,WD_1_1_1,P,
%SW_IN,SW_OUT,LW_IN,LW_OUT,T_CANOPY_1_1_1,T_CANOPY_2_1_1,PPFD_IN_PI_F,PPFD_OUT,PPFD_IN,
%SWC_PI_1_1_A,SWC_PI_1_2_A,SWC_PI_1_3_A,SWC_PI_1_4_A,SWC_PI_1_5_A,SWC_PI_1_6_A,SWC_PI_1_7_A,SWC_PI_1_8_A,TS_PI_1_1_A,
%TS_PI_1_2_A,TS_PI_1_3_A,TS_1_4_1,TS_PI_1_5_A,TS_PI_1_6_A,TS_PI_1_7_A,TS_PI_1_8_A,TS_PI_1_4_A


%%%%%%%%%%%  Defining Extracting Variables
Var_to_extract={'TA';'P';'WD';'WS';'NEE_PI';'FC';'H';'LE';'RH';'PA';'G';'VPD_PI';...
    'CO2_PI_2';'CO2_2';...
    'TS_PI_1';'TS_PI_2';'TS_1';'TS_2';...
    'SWC_PI_1';'SWC_PI_2';'SWC_1';'SWC_2';...
    'NETRAD';'PPFD_IN';'SW_IN';'SW_DIF';'PPFD_OUT';'SW_OUT';'LW_IN';'LW_OUT';...
    'RECO_PI';'APAR';'FAPAR';'PPFD_DIF';'GPP_PI'};

for k=1:length(Var_to_extract)
    II=strcmp(Var_List,Var_to_extract{k});
    pos = find(II==1);
    if isempty(pos)
        PV(k)=NaN;
    else
        PV(k)=pos;
    end
end
clear II k pos

%%% If you want to reduce the extent of the time series 
x1=1;
x2=262990; 
DATA=DATA(x1:x2,:); 

d1=datenum(num2str(DATA(1,2)),'yyyymmddHHMM');
d2=datenum(num2str(DATA(end,2)),'yyyymmddHHMM');
% Hourly Date
Date = d1:1/24:d2; Date=Date';



%%%%%%  Convert data: half-hourly to hourly
dt = 0.5;
n=length(DATA(:,1)); fr=1/dt;  m=floor(n/fr);
j=0; DATA1h=NaN*ones(m,length(Var_to_extract));
%
for j=1:length(Var_to_extract)
    if not(isnan(PV(j)))
        V=DATA(:,PV(j));
        V=reshape(V(1:m*fr),fr,m);
        if(j==2) %% Precipitation
            X= nansum(V);
        else
            X= nanmean(V);
        end
        DATA1h(:,j)=X;
        clear V X
    end
end
%
DATA=DATA1h;

%%%%%%%

Ta =DATA(:,1);    Ta(Ta<-60)=NaN; Ta(Ta>70)=NaN;                            %
Pr =DATA(:,2);    Pr(Pr<0)=NaN; Pr(Pr>220)=NaN;                            %
Wdir = DATA(:,3);
Ws = DATA(:,4);  Ws(Ws<0)=NaN; Ws(Ws>100)=NaN;                             %
NEE = DATA(:,5); NEE(NEE<-500)=NaN; NEE(NEE>500)=NaN;                      %
FC = DATA(:,6);
H =DATA(:,7);   H(H<-1550)=NaN; H(H>1500)=NaN;
LE = DATA(:,8);  LE(LE<-250)=NaN; LE(LE>1500)=NaN;
U = DATA(:,9);   U(U<0)=NaN; U(U>100)=NaN;                                 % 
Pre = DATA(:,10)*10; Pre(Pre<0)=NaN;                                        %
G = DATA(:,11);  G(G<-250)=NaN; G(G>1500)=NaN;
VPD = DATA(:,12);  VPD(VPD<0)=NaN; VPD(VPD>70)=NaN;
CO2 = DATA(:,14); CO2(CO2<0)=NaN; CO2(CO2>900)=NaN;
if nansum(CO2)==0
    CO2 = DATA(:,13); CO2(CO2<0)=NaN; CO2(CO2>900)=NaN;
end
Tsoil1= DATA(:,17) ; Tsoil1(Tsoil1<-50)=NaN; Tsoil1(Tsoil1>90)=NaN;         %
Tsoil2= DATA(:,18) ; Tsoil2(Tsoil2<-50)=NaN; Tsoil2(Tsoil2>90)=NaN;         %
if nansum(Tsoil1)==0
    Tsoil1= DATA(:,15) ; Tsoil1(Tsoil1<-50)=NaN; Tsoil1(Tsoil1>90)=NaN;         %
    Tsoil2= DATA(:,16) ; Tsoil2(Tsoil2<-50)=NaN; Tsoil2(Tsoil2>90)=NaN;         %
end
SWC1 = DATA(:,21);  SWC1(SWC1<0)=NaN; SWC1(SWC1>100)=NaN;                   %
SWC2 = DATA(:,22);  SWC2(SWC2<0)=NaN; SWC2(SWC2>100)=NaN;                   %
if nansum(SWC1)==0
    SWC1 = DATA(:,19);  SWC1(SWC1<0)=NaN; SWC1(SWC1>100)=NaN;                   %
    SWC2 = DATA(:,20);  SWC2(SWC2<0)=NaN; SWC2(SWC2>100)=NaN;                   %
end
Rn =DATA(:,23); Rn(Rn<-1550)=NaN; Rn(Rn>1000)=NaN;                          % 
PPFD = DATA(:,24); PPFD(PPFD<0)=NaN; PPFD(PPFD>3000)=NaN;
Rsw = DATA(:,25);  Rsw(Rsw<-10)=NaN; Rsw(Rsw>2000)=NaN; Rsw(Rsw<0)=0;
Rdif = DATA(:,26);  Rdif(Rdif<-10)=NaN; Rdif(Rdif>2000)=NaN; Rdif(Rdif<0)=0;
PPFD_out = DATA(:,27); PPFD_out(PPFD_out<0)=NaN; PPFD_out(PPFD_out>3000)=NaN;
Rsw_out = DATA(:,28);  Rsw_out(Rsw_out<-10)=NaN; Rsw_out(Rsw_out>2000)=NaN; Rsw_out(Rsw_out<0)=0;
Latm = DATA(:,29); Latm(Latm<-600)=NaN; Latm(Latm>600)=NaN;
Lup = DATA(:,30); Lup(Lup<-600)=NaN; Lup(Lup>600)=NaN;
Reco = DATA(:,31); Reco(Reco<-10)=NaN; Reco(Reco>1500)=NaN;
APAR = DATA(:,32);
fPAR = DATA(:,33);
PPFD_dif = DATA(:,34);
GPP = DATA(:,35);  GPP(GPP<-10)=NaN;
if nansum(GPP)==0
    GPP = - NEE + Reco;
end
%%%
U=U/100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Figure (pre smoothing) ---------------------------------

figure(101)
subplot(4,2,1)
plot(Date,Ta); hold on
datetick('x',11); title('Ta')
subplot(4,2,2)
plot(Date,Pr); hold on
datetick('x',11); title(strcat('Pr:',num2str(nanmean(Pr)*8760)))
subplot(4,2,3)
plot(Date,Rsw); hold on
datetick('x',11); title(strcat('Rsw:',num2str(nanmean(Rsw)*1)))
subplot(4,2,4)
plot(Date,U);hold on
datetick('x',11); title('U')
subplot(4,2,5)
plot(Date,Pre);hold on
datetick('x',11); title('Pre')
subplot(4,2,6)
plot(Date,Ws);hold on
datetick('x',11); title('Ws')
subplot(4,2,7)
plot(Date,CO2);hold on
datetick('x',11); title('Ca')
subplot(4,2,8)
plot(Date,Latm);hold on
datetick('x',11); title('Lwd')


figure(102)
subplot(3,2,1)
plot(Date,H); hold on
datetick('x',11); title('H')
subplot(3,2,2)
plot(Date,LE); hold on
datetick('x',11); title('\lambda E')
subplot(3,2,3)
plot(Date,G); hold on
datetick('x',11); title('G')
subplot(3,2,4)
plot(Date,GPP);hold on
datetick('x',11); title('GPP')
subplot(3,2,5)
plot(Date,NEE,'r');hold on
datetick('x',11); title('NEE')
subplot(3,2,6)
plot(Date,Rn);hold on
datetick('x',11); title('Rn')

figure(103)
subplot(3,2,1)
plot(Date,Tsoil1); hold on
hold on
plot(Date,Tsoil2,'--'); hold on
datetick('x',11); title('Ts')
subplot(3,2,2)
plot(Date,SWC1); hold on
plot(Date,SWC2,'--');
datetick('x',11); title('SWC')
subplot(3,2,3)
plot(Date,Reco);hold on
datetick('x',11); title('R_{eco}')
subplot(3,2,4)
plot(Date,PPFD);hold on
datetick('x',11); title('PPFD')
subplot(3,2,5)
plot(Date,Lup);hold on
datetick('x',11); title('Lw_{up}')
subplot(3,2,6)
plot(Date,Rsw_out);hold on
datetick('x',11);  title('Rsw_{up}')


r=0;
Yrs=year(Date);
for i=min(Yrs):max(Yrs)
    r=r+1;
    Pr_ann=Pr(find(Yrs==i));
    Pr_yr(r)= nansum(Pr_ann);
    Pr_yr2(r)= nanmean(Pr_ann)*8766;
    %%%%%
end
figure(104)
plot(min(Yrs):max(Yrs),Pr_yr,'-ob');
hold on
plot(min(Yrs):max(Yrs),Pr_yr2,'--db');
legend('nansum(Pr_{ann})','nanmean(Pr_{ann})*8766')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% FILL GAPS -- INTERPOLATION NaNs
Pr(isnan(Pr))=0;

Pre(isnan(Pre))=nanmean(Pre);
if nansum(Pre)==0
    p=1013.25*exp(-Zbas/8434.5); %%%
    Pre2= p*ones(1,length(Date)); clear p
    Pre(isnan(Pre))=Pre2(isnan(Pre));
end

dT =[0; diff(Ta)];
Ta(abs(dT)>10.0)=NaN;
dT =[0; diff(Ta)];
Ta(isnan(dT))=NaN;
dT =[0; diff(Ta)];
Ta(isnan(dT))=NaN;
dWs =[0; diff(Ws)];
Ws(abs(dWs)>9)=NaN;
dWs =[0; diff(Ws)];
Ws(isnan(dWs))=NaN;
dWs =[0; diff(Ws)];
Ws(isnan(dWs))=NaN;
dCO2 =[0; diff(CO2)];
CO2(abs(dCO2)>40)=NaN;
dCO2 =[0; diff(CO2)];
CO2(isnan(dCO2))=NaN;
dCO2 =[0; diff(CO2)];
CO2(isnan(dCO2))=NaN;

% SMOOTHING/INTERPOLATION
%---------------------------------------------
cd('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data')
[Ta]=Interp_Smooth_Cycle(Ta,Date);
[Ws]=Interp_Smooth_Cycle(Ws,Date);
[U]=Interp_Smooth_Cycle(U,Date);
[Pr]=Interp_Smooth_Cycle_Pr(Pr,Date);
[CO2]=Interp_Smooth_Cycle(CO2,Date);
[Latm]=Interp_Smooth_Cycle(Latm,Date);
cd(cur_dir);

U(U<0.04)=NaN;
U(isnan(U))=interp1(Date(not(isnan(U))),U(not(isnan(U))),Date(isnan(U)));
U(isnan(U))=nanmean(U);

Ta(isnan(Ta))=interp1(Date(not(isnan(Ta))),Ta(not(isnan(Ta))),Date(isnan(Ta)));
Ta(isnan(Ta))=nanmean(Ta);
%%%%%%%%
Ws(isnan(Ws))=interp1(Date(not(isnan(Ws))),Ws(not(isnan(Ws))),Date(isnan(Ws)));
Ws(isnan(Ws))=nanmean(Ws);

Latm(isnan(Latm))=interp1(Date(not(isnan(Latm))),Latm(not(isnan(Latm))),Date(isnan(Latm)));
Latm(isnan(Latm))=nanmean(Latm);

%CO2(isnan(CO2))=nanmean(CO2);

esat=611*exp(17.27*Ta./(237.3+Ta)); %% [Pa] Vapor pressure saturation
ea=esat.*U; %% Vapor Pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=17.27; b=237.7;
xr= a*Ta./(b+Ta) + log(U);
Tdew =b*xr./(a-xr);
clear a b xr

% Figure (post interpolation) ---------------------------------------------
%

figure(101)
subplot(4,2,1)
plot(Date,Ta); hold on
datetick('x',11); title('Ta')
subplot(4,2,2)
plot(Date,Pr); hold on
datetick('x',11); title(strcat('Pr:',num2str(nanmean(Pr)*8760)))
subplot(4,2,3)
plot(Date,Rsw); hold on
datetick('x',11); title(strcat('Rsw:',num2str(nanmean(Rsw)*1)))
subplot(4,2,4)
plot(Date,U);hold on
datetick('x',11); title('U')
subplot(4,2,5)
plot(Date,Pre);hold on
datetick('x',11); title('Pre')
subplot(4,2,6)
plot(Date,Ws);hold on
datetick('x',11); title('Ws')
subplot(4,2,7)
plot(Date,CO2);hold on
datetick('x',11); title('Ca')
subplot(4,2,8)
plot(Date,Latm);hold on
datetick('x',11); title('Lwd')


figure(102)
subplot(3,2,1)
plot(Date,H); hold on
datetick('x',11); title('H')
subplot(3,2,2)
plot(Date,LE); hold on
datetick('x',11); title('\lambda E')
subplot(3,2,3)
plot(Date,G); hold on
datetick('x',11); title('G')
subplot(3,2,4)
plot(Date,GPP);hold on
datetick('x',11); title('GPP')
subplot(3,2,5)
plot(Date,NEE,'r');hold on
datetick('x',11); title('NEE')
subplot(3,2,6)
plot(Date,Rn);hold on
datetick('x',11); title('Rn')

figure(103)
subplot(3,2,1)
plot(Date,Tsoil1); hold on
hold on
plot(Date,Tsoil2,'--'); hold on
datetick('x',11); title('Ts')
subplot(3,2,2)
plot(Date,SWC1); hold on
plot(Date,SWC2,'--');
datetick('x',11); title('SWC')
subplot(3,2,3)
plot(Date,Reco);hold on
datetick('x',11); title('R_{eco}')
subplot(3,2,4)
plot(Date,PPFD);hold on
datetick('x',11); title('PPFD')
subplot(3,2,5)
plot(Date,Lup);hold on
datetick('x',11); title('Lw_{up}')
subplot(3,2,6)
plot(Date,Rsw_out);hold on
datetick('x',11);  title('Rsw_{up}')


r=0;
Yrs=year(Date);
for i=min(Yrs):max(Yrs)
    r=r+1;
    Pr_ann=Pr(find(Yrs==i));
    Pr_yr(r)= nansum(Pr_ann);
    Pr_yr2(r)= nanmean(Pr_ann)*8766;
    %%%%%
end
figure(104)
plot(min(Yrs):max(Yrs),Pr_yr,'-ob');
hold on
plot(min(Yrs):max(Yrs),Pr_yr2,'--db');
legend('nansum(Pr_{ann})','nanmean(Pr_{ann})*8766')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%% 
save(['Res_',Location_Name,'.mat'],'Date','DeltaGMT','Lon','Lat','Zbas','Ta','Rsw','esat','ea','Tdew',...
    'PPFD','Latm','VPD','Pre','Pr','Ws','U','CO2','Tsoil1','Tsoil2','Tsoil3','Tsoil4','Tsoil5','Tsoil6','SWC1','SWC2','SWC3','SWC4','SWC5','SWC6','SWC7',... 
'Rn','LE','LE_CORR','H','H_CORR','G','NEE1','NEE2','NEE3','NEE4','Reco1','Reco2','Reco3','Reco4','GPP1','GPP2','GPP3','GPP4');  

clearvars -except Location_Name cur_dir


[SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,Rsws,t_bef,t_aft]=Automatic_Radiation_Partition(cur_dir,Location_Name,1); 


load(['Res_',Location_Name,'.mat']);

%%% Sun Integration [0.0 1.0]


save(['Data_',Location_Name,'_run.mat'],'Date','Pr','Ta','Ws','U','ea','SAD1','SAD2','SAB1','SAB2','N','Latm','Tdew','esat','PARB','PARD','Pre',...
   'DeltaGMT','Lon','Lat','Zbas','CO2','t_bef','t_aft');
