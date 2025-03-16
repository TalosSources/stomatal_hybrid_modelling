d%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Run Section by Section!!!
%--------------------------------
%The following code is reading FLUXNET2015 half hourly (HH) data %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Site description 	
%%%		 IT-Ro2: Roccarespampani 2 // Dec. Forest 
Lon          = 11.9209;  
Lat          = 42.3903		;
DeltaGMT     = 1;
Zbas         = 160  ; 
Location_Name ='Roccarespampani2'; 
		
		
cur_dir = cd;
file_dir = ('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_IT-Ro2_FLUXNET2015_SUBSET_2002-2012_1-3');

file_name = 'FLX_IT-Ro2_FLUXNET2015_SUBSET_HH_2002-2012_1-3.csv';

%%
A = importdata(strcat(file_dir,'\',file_name),',',1);
Var_List = A.textdata; DATA=A.data;  clear A
DATA(DATA==-9999)=NaN;

%%% 227 fields 
%%% TIMESTAMP_START,TIMESTAMP_END,TA_F_MDS,TA_F_MDS_QC,TA_ERA,TA_F,TA_F_QC,SW_IN_POT,SW_IN_F_MDS,SW_IN_F_MDS_QC,SW_IN_ERA,SW_IN_F,SW_IN_F_QC,LW_IN_F_MDS,LW_IN_F_MDS_QC,LW_IN_ERA,
%%%LW_IN_F,LW_IN_F_QC,LW_IN_JSB,LW_IN_JSB_QC,LW_IN_JSB_ERA,LW_IN_JSB_F,LW_IN_JSB_F_QC,VPD_F_MDS,VPD_F_MDS_QC,VPD_ERA,VPD_F,VPD_F_QC,PA,PA_ERA,PA_F,
%%%PA_F_QC,P,P_ERA,P_F,P_F_QC,WS,WS_ERA,WS_F,WS_F_QC,WD,USTAR,RH,NETRAD,PPFD_IN,SW_OUT,CO2_F_MDS,CO2_F_MDS_QC,TS_F_MDS_1,TS_F_MDS_2,TS_F_MDS_3,
%%TS_F_MDS_1_QC,TS_F_MDS_2_QC,TS_F_MDS_3_QC,SWC_F_MDS_1,SWC_F_MDS_1_QC,G_F_MDS,G_F_MDS_QC,LE_F_MDS,LE_F_MDS_QC,LE_CORR,LE_CORR_25,LE_CORR_75,
%%LE_RANDUNC,LE_RANDUNC_METHOD,LE_RANDUNC_N,LE_CORR_JOINTUNC,H_F_MDS,H_F_MDS_QC,H_CORR,H_CORR_25,H_CORR_75,H_RANDUNC,H_RANDUNC_METHOD,H_RANDUNC_N,
%%H_CORR_JOINTUNC,EBC_CF_N,EBC_CF_METHOD,NIGHT,NEE_CUT_REF,NEE_VUT_REF,NEE_CUT_REF_QC,NEE_VUT_REF_QC,NEE_CUT_REF_RANDUNC,NEE_VUT_REF_RANDUNC,
%%NEE_CUT_REF_RANDUNC_METHOD,NEE_VUT_REF_RANDUNC_METHOD,NEE_CUT_REF_RANDUNC_N,NEE_VUT_REF_RANDUNC_N,NEE_CUT_REF_JOINTUNC,NEE_VUT_REF_JOINTUNC,NEE_CUT_USTAR50,
%%NEE_VUT_USTAR50,NEE_CUT_USTAR50_QC,NEE_VUT_USTAR50_QC,NEE_CUT_USTAR50_RANDUNC,NEE_VUT_USTAR50_RANDUNC,NEE_CUT_USTAR50_RANDUNC_METHOD,
%NEE_VUT_USTAR50_RANDUNC_METHOD,NEE_CUT_USTAR50_RANDUNC_N,NEE_VUT_USTAR50_RANDUNC_N,NEE_CUT_USTAR50_JOINTUNC,NEE_VUT_USTAR50_JOINTUNC,NEE_CUT_MEAN,NEE_VUT_MEAN,
%NEE_CUT_MEAN_QC,NEE_VUT_MEAN_QC,NEE_CUT_SE,NEE_VUT_SE,NEE_CUT_05,NEE_CUT_16,NEE_CUT_25,NEE_CUT_50,NEE_CUT_75,NEE_CUT_84,NEE_CUT_95,NEE_VUT_05,
%NEE_VUT_16,NEE_VUT_25,NEE_VUT_50,NEE_VUT_75,NEE_VUT_84,NEE_VUT_95,NEE_CUT_05_QC,NEE_CUT_16_QC,NEE_CUT_25_QC,NEE_CUT_50_QC,NEE_CUT_75_QC,
%NEE_CUT_84_QC,NEE_CUT_95_QC,NEE_VUT_05_QC,NEE_VUT_16_QC,NEE_VUT_25_QC,NEE_VUT_50_QC,NEE_VUT_75_QC,NEE_VUT_84_QC,NEE_VUT_95_QC,
%RECO_NT_VUT_REF,RECO_NT_VUT_USTAR50,RECO_NT_VUT_MEAN,RECO_NT_VUT_SE,RECO_NT_VUT_05,RECO_NT_VUT_16,RECO_NT_VUT_25,RECO_NT_VUT_50,
%RECO_NT_VUT_75,RECO_NT_VUT_84,RECO_NT_VUT_95,RECO_NT_CUT_REF,RECO_NT_CUT_USTAR50,RECO_NT_CUT_MEAN,RECO_NT_CUT_SE,RECO_NT_CUT_05,
%RECO_NT_CUT_16,RECO_NT_CUT_25,RECO_NT_CUT_50,RECO_NT_CUT_75,RECO_NT_CUT_84,RECO_NT_CUT_95,GPP_NT_VUT_REF,GPP_NT_VUT_USTAR50,
%GPP_NT_VUT_MEAN,GPP_NT_VUT_SE,GPP_NT_VUT_05,GPP_NT_VUT_16,GPP_NT_VUT_25,GPP_NT_VUT_50,GPP_NT_VUT_75,GPP_NT_VUT_84,GPP_NT_VUT_95,
%GPP_NT_CUT_REF,GPP_NT_CUT_USTAR50,GPP_NT_CUT_MEAN,GPP_NT_CUT_SE,GPP_NT_CUT_05,GPP_NT_CUT_16,GPP_NT_CUT_25,GPP_NT_CUT_50,
%GPP_NT_CUT_75,GPP_NT_CUT_84,GPP_NT_CUT_95,RECO_DT_VUT_REF,RECO_DT_VUT_USTAR50,RECO_DT_VUT_MEAN,RECO_DT_VUT_SE,RECO_DT_VUT_05,
%RECO_DT_VUT_16,RECO_DT_VUT_25,RECO_DT_VUT_50,RECO_DT_VUT_75,RECO_DT_VUT_84,RECO_DT_VUT_95,RECO_DT_CUT_REF,RECO_DT_CUT_USTAR50,
%RECO_DT_CUT_MEAN,RECO_DT_CUT_SE,RECO_DT_CUT_05,RECO_DT_CUT_16,RECO_DT_CUT_25,RECO_DT_CUT_50,RECO_DT_CUT_75,RECO_DT_CUT_84,
%RECO_DT_CUT_95,GPP_DT_VUT_REF,GPP_DT_VUT_USTAR50,GPP_DT_VUT_MEAN,GPP_DT_VUT_SE,GPP_DT_VUT_05,GPP_DT_VUT_16,GPP_DT_VUT_25,
%GPP_DT_VUT_50,GPP_DT_VUT_75,GPP_DT_VUT_84,GPP_DT_VUT_95,GPP_DT_CUT_REF,GPP_DT_CUT_USTAR50,GPP_DT_CUT_MEAN,GPP_DT_CUT_SE,
%GPP_DT_CUT_05,GPP_DT_CUT_16,GPP_DT_CUT_25,GPP_DT_CUT_50,GPP_DT_CUT_75,GPP_DT_CUT_84,GPP_DT_CUT_95,RECO_SR,RECO_SR_N

%%%%%%%%%%%  Defining Extracting Variables 
Var_to_extract={'TA_F';'SW_IN_F';'PPFD_IN';'LW_IN_F';'VPD_F';'PA_F';'P_F';'WS_F';'RH';'CO2_F_MDS';...
    'TS_F_MDS_1';'TS_F_MDS_2';'TS_F_MDS_3';'TS_F_MDS_4';'TS_F_MDS_5';'TS_F_MDS_6';...
    'SWC_F_MDS_1';'SWC_F_MDS_2';'SWC_F_MDS_3';'SWC_F_MDS_4';'SWC_F_MDS_5';'SWC_F_MDS_6';'SWC_F_MDS_7';...
    'NETRAD';'LE_F_MDS';'LE_F_MDS_QC';'LE_CORR';'H_F_MDS';'H_F_MDS_QC';'H_CORR';'G_F_MDS'; 'NEE_CUT_REF';'NEE_VUT_REF';...
    'NEE_CUT_50';'NEE_VUT_50';'RECO_NT_VUT_REF';'RECO_NT_VUT_50 ';'RECO_DT_VUT_REF';'RECO_DT_VUT_50';...
    'GPP_NT_VUT_REF';'GPP_NT_VUT_50';'GPP_DT_VUT_REF';'GPP_DT_VUT_50'};

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
        if(j==7) %% Precipitation 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th1=3; %% Delta W/m2 / h  
th2=0.2; %% Delta uCO2/m2 s / h 
%%% Cleaning of repetead Variables badly interpolated 
for i=[24 25 28 31]
        V= DATA(:,i);
        V3=V;
        for j=336:168:length(V)
            if sum(abs(V(j-335:j-168)-V(j-167:j)))< th1*168
                V3(j-167:j)=NaN;
            end
        end
        V=V3; V2=V;
        for j=144:72:length(V)
            if sum(abs(V(j-143:j-72)-V(j-71:j)))< th1*72
                V2(j-71:j)=NaN;
            end
        end
        V=V2; V1=V;
        for j=48:24:length(V)
            if sum(abs(V(j-47:j-24)-V(j-23:j)))< th1*24
                V2(j-23:j)=NaN;
            end
        end
        DATA(:,i)=V1;
        clear V2 V V3 V1 
end
for i=32:43
        V= DATA(:,i);
        V3=V;
        for j=336:168:length(V)
            if sum(abs(V(j-335:j-168)-V(j-167:j)))< th2*168
                V3(j-167:j)=NaN;
            end
        end
        V=V3; V2=V;
        for j=144:72:length(V)
            if sum(abs(V(j-143:j-72)-V(j-71:j)))< th2*72
                V2(j-71:j)=NaN;
            end
        end
        V=V2; V1=V;
        for j=48:24:length(V)
            if sum(abs(V(j-47:j-24)-V(j-23:j)))< th2*24
                V2(j-23:j)=NaN;
            end
        end
        DATA(:,i)=V1;
        clear V2 V V3 V1 
end    

%%%%%%%%
%Ta, Rsw PPFD Latm VPD  Pre  Pr Ws RH CO2 Tsoil1 Tsoil2 Tsoil3 Tsoil4 Tsoil5 Tsoil6 SWC1 SWC2 SWC3 SWC4 SWC5 SWC6 SWC7 
%Rn LE LE_QC LE_CORR H H_QC H_CORR G NEE1 NEE2 NEE3 NEE4  Reco1 Reco2 Reco3 Reco4 GPP1 GPP2 GPP3 GPP4 

Ta =DATA(:,1);    Ta(Ta<-60)=NaN; Ta(Ta>70)=NaN;                            % TA_F - degC
Rsw = DATA(:,2);  Rsw(Rsw<-10)=NaN; Rsw(Rsw>2000)=NaN; Rsw(Rsw<0)=0;       % SW_IN_F - W/m2 - Incoming global solar radiation (W m-2) 
PPFD = DATA(:,3); PPFD(PPFD<-10)=NaN; PPFD(PPFD>3000)=NaN;   PPFD(PPFD<0)=0;    % PPFD_IN - W/m2 
Latm = DATA(:,4); Latm(Latm<-600)=NaN; Latm(Latm>600)=NaN; 
VPD = DATA(:,5);  VPD(VPD<0)=NaN; VPD(VPD>70)=NaN;                         % VPD_F - hPa
Pre = DATA(:,6)*10; Pre(Pre<0)=NaN;                                        % PA_F - kPa -> mbar
Pr =DATA(:,7);    Pr(Pr<0)=NaN; Pr(Pr>220)=NaN;                            % P_F - mm (35)
Ws = DATA(:,8);  Ws(Ws<0)=NaN; Ws(Ws>100)=NaN;                             % WS_F - m/s
U = DATA(:,9);   U(U<0)=NaN; U(U>100)=NaN;                                 % RH - %
CO2 = DATA(:,10); CO2(CO2<0)=NaN; CO2(CO2>900)=NaN;                         % CO2 concentration (umol mol-1) 
Tsoil1= DATA(:,11) ; Tsoil1(Tsoil1<-50)=NaN; Tsoil1(Tsoil1>90)=NaN;         % TS_F_MDS_1 - degC - Soil temperature 1 (degrees C) 
Tsoil2= DATA(:,12) ; Tsoil2(Tsoil2<-50)=NaN; Tsoil2(Tsoil2>90)=NaN;         % TS_F_MDS_2 - degC - Soil temperature  2 (degrees C) 
Tsoil3= DATA(:,13) ; Tsoil3(Tsoil3<-50)=NaN; Tsoil3(Tsoil3>90)=NaN;         % TS_F_MDS_3 - degC - Soil temperature 3 (degrees C) 
Tsoil4= DATA(:,14) ; Tsoil4(Tsoil4<-50)=NaN; Tsoil4(Tsoil4>90)=NaN;         % TS_F_MDS_4 - degC - Soil temperature 4 (degrees C)
Tsoil5= DATA(:,15) ; Tsoil5(Tsoil5<-50)=NaN; Tsoil5(Tsoil5>90)=NaN;         % TS_F_MDS_5 - degC - Soil temperature 4 (degrees C)
Tsoil6= DATA(:,16) ; Tsoil6(Tsoil6<-50)=NaN; Tsoil6(Tsoil6>90)=NaN;         % TS_F_MDS_6 - degC - Soil temperature 4 (degrees C)
SWC1 = DATA(:,17);  SWC1(SWC1<0)=NaN; SWC1(SWC1>100)=NaN;                   % SWC_F_MDS_1 - % - Soil moisture content at depth 1
SWC2 = DATA(:,18);  SWC2(SWC2<0)=NaN; SWC2(SWC2>100)=NaN;                   % SWC_F_MDS_2 - % - Soil moisture content at depth 2
SWC3 = DATA(:,19);  SWC3(SWC3<0)=NaN; SWC3(SWC3>100)=NaN;                   % SWC_F_MDS_3 - % - Soil moisture content at depth 3
SWC4 = DATA(:,20);  SWC4(SWC4<0)=NaN; SWC4(SWC4>100)=NaN;                   % SWC_F_MDS_4 - % - Soil moisture content at depth 4
SWC5 = DATA(:,21);  SWC5(SWC5<0)=NaN; SWC5(SWC5>100)=NaN;                   % SWC_F_MDS_5 - % - Soil moisture content at depth 5
SWC6 = DATA(:,22);  SWC6(SWC6<0)=NaN; SWC6(SWC6>100)=NaN;                   % SWC_F_MDS_6 - % - Soil moisture content at depth 6
SWC7 = DATA(:,23);  SWC7(SWC7<0)=NaN; SWC7(SWC7>100)=NaN;                   % SWC_F_MDS_7 - % - Soil moisture content at depth 7
%%%%
Rn =DATA(:,24); Rn(Rn<-1550)=NaN; Rn(Rn>1000)=NaN;                          % NETRAD - W/m2
LE = DATA(:,25);  LE(LE<-250)=NaN; LE(LE>1500)=NaN;                         % LE_F_MDS - W/m2
LE_QC=DATA(:,26);%0 = measured; 1 = good quality gapfill; 2 = medium; 3 = poor
LE_CORR = DATA(:,27); %%% LE corrected for energy balance closure
H =DATA(:,28);   H(H<-1550)=NaN; H(H>1500)=NaN;                             % H_F_MDS - W/m2
H_QC=DATA(:,29);%0 = measured; 1 = good quality gapfill; 2 = medium; 3 = poor
H_CORR = DATA(:,30); %%% H corrected for energy balance closure
G = DATA(:,31);  G(G<-250)=NaN; G(G>1500)=NaN;                              % G_F_MDS - W/m2 - Soil heat flux (W m-2) 
NEE1 = DATA(:,32); NEE1(NEE1<-500)=NaN; NEE1(NEE1>500)=NaN;                      % NEE_CUT_REF - umol/m2/s - NEE of CO2 (with storage) (umol m-2 s-1)
NEE2 = DATA(:,33); NEE2(NEE2<-500)=NaN; NEE2(NEE2>500)=NaN;                      % NEE_VUT_REF - umol/m2/s - NEE of CO2 (with storage) (umol m-2 s-1)
NEE3 = DATA(:,34); NEE3(NEE3<-500)=NaN; NEE3(NEE3>500)=NaN;                      
NEE4 = DATA(:,35); NEE4(NEE4<-500)=NaN; NEE4(NEE4>500)=NaN;                   
Reco1 = DATA(:,36); Reco1(Reco1<-10)=NaN; Reco1(Reco1>1500)=NaN;  
Reco2 =DATA(:,37);  Reco2(Reco2<-10)=NaN; Reco2(Reco2>1500)=NaN;  
Reco3 =DATA(:,38);  Reco3(Reco3<-10)=NaN; Reco3(Reco3>1500)=NaN;  
Reco4 =DATA(:,39);  Reco4(Reco4<-10)=NaN; Reco4(Reco4>1500)=NaN;  
GPP1 = DATA(:,40);  GPP1(GPP1<-10)=NaN;                                         % GPP_NT_VUT_REF (nighttime) - umo/m2/s - Gross primary production (umol/m2/s)
GPP2 = DATA(:,41);  GPP2(GPP2<-10)=NaN;                                      % GPP_DT_VUT_REF (daytime) - umo/m2/s - Gross primary production (umol/m2/s)
GPP3 = DATA(:,42);  GPP3(GPP3<-10)=NaN;   
GPP4 = DATA(:,43);  GPP4(GPP4<-10)=NaN;   

U=U/100;
%esat=611*exp(17.27*Ta./(237.3+Ta)); %% [Pa] Vapor pressure saturation
%ea = esat- VPD*100;   U2 = ea./esat; U(isnan(U))=U2(isnan(U)); 
LE(LE_QC==3)=NaN; H(H_QC==3)=NaN; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
plot(Date,GPP1);hold on
plot(Date,GPP2,'--')
plot(Date,GPP3,'--')
plot(Date,GPP4,'--')
datetick('x',11); title('GPP')
subplot(3,2,5)
plot(Date,NEE1,'r');hold on
plot(Date,NEE2,'--r');hold on
plot(Date,NEE3,'--r');hold on
plot(Date,NEE4,'--r');hold on
datetick('x',11); title('NEE')
subplot(3,2,6)
plot(Date,Rn);hold on
datetick('x',11); title('Rn')

figure(103)
subplot(2,2,1)
plot(Date,Tsoil1); hold on
hold on 
plot(Date,Tsoil2); hold on
plot(Date,Tsoil3); hold on
datetick('x',11); title('Ts')
subplot(2,2,2)
plot(Date,SWC1); hold on
plot(Date,SWC2);
plot(Date,SWC3);
%plot(Date,SWC4);
datetick('x',11); title('SWC')
subplot(2,2,3)
plot(Date,Reco1);hold on
plot(Date,Reco2,'--');hold on
plot(Date,Reco3,'--');hold on
plot(Date,Reco4,'--');hold on
datetick('x',11); title('R_{eco}')
subplot(2,2,4)
plot(Date,PPFD);hold on
datetick('x',11); title('PPFD')

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

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% FILL GAPS -- INTERPOLATION NaNs
%Pr(isnan(Pr))=0;

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
datetick('x',11); title('Latm')


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
plot(min(Yrs):max(Yrs),Pr_yr,'-or');
hold on
plot(min(Yrs):max(Yrs),Pr_yr2,'--dr');
legend('nansum(Pr_{ann})','nanmean(Pr_{ann})*8766')
%--------------------------------------------------------------------------
%%-
%%
%%%%%%%% 
save(['Res_',Location_Name,'.mat'],'Date','DeltaGMT','Lon','Lat','Zbas','Ta','Rsw','esat','ea','Tdew',...
    'PPFD','Latm','VPD','Pre','Pr','Ws','U','CO2','Tsoil1','Tsoil2','Tsoil3','Tsoil4','Tsoil5','Tsoil6','SWC1','SWC2','SWC3','SWC4','SWC5','SWC6','SWC7',... 
'Rn','LE','LE_CORR','H','H_CORR','G','NEE1','NEE2','NEE3','NEE4','Reco1','Reco2','Reco3','Reco4','GPP1','GPP2','GPP3','GPP4');  

clearvars -except Location_Name cur_dir

cd('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data')
[SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,Rsws,t_bef,t_aft]=Automatic_Radiation_Partition(cur_dir,Location_Name,1); 
cd(cur_dir)

%cd('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data')
%[SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,Rsws,t_bef,t_aft]=Automatic_Radiation_Partition_I(Date,Lat,Lon,Zbas,DeltaGMT,Pr,Tdew,Rsw,1);
%cd(cur_dir);

load(['Res_',Location_Name,'.mat']);

%%% Sun Integration [0.0 1.0]


save(['Data_',Location_Name,'_run.mat'],'Date','Pr','Ta','Ws','U','ea','SAD1','SAD2','SAB1','SAB2','N','Latm','Tdew','esat','PARB','PARD','Pre',...
   'DeltaGMT','Lon','Lat','Zbas','CO2','t_bef','t_aft');

