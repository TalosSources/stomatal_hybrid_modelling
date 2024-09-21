%%%%%%%%
%%% dB/dt = -B*tt + ANPPw
%%% B= ANPPw/tt ; 
%% 1/tt =  B/ANPPw; 
a = 120:1:1500;%% ANPP [gDM/m2 yr] Aboveground NPP 
aw = a/3;  %% aboveground wood allocation 
%Es_AGBfit =  [Mg DM / ha] Expected Maximum Standing Biomass 
Es_AGBfit = (30.34.*a.^0.3735-178.3); %%%  Fit of 80% envelope 
%Es_AGBfit = (15.42 +22.96*a -0.498*a.^2); %
%%%%
tt = 100*Es_AGBfit./aw; %%%  [yr] 
%%%%%%%
figure(101)
plot(a,tt,'o'); 
xlabel('ANPP [gDM/m2 yr]'); ylabel('Wood turnover rate [yr]')




%%%%%%%%%%%%%%%%%%%%
a=433*2; %% gDM m-2 yr-1 
aw=204*0.74*2;  %% gDM m-2 yr-1 
tt=66;  %% yr
AWOOD=0; 
for yy=2:500
    AWOOD(yy)=AWOOD(yy-1)*(1-1/tt)+aw; 
end 

figure(102)
plot(1:500,AWOOD/100)
xlabel('Year [yr]'); ylabel('Aboveground Biomass [MgDM ha-1]')