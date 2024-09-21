%%%%%%%%%%%%%%%%
clear
OPT_CASE =105;
cur_dir2=cd;
switch OPT_CASE
    case 1
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Renon\Res_Renon.mat')
    case 2
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Valles_Caldera_Mix_Con\Res_VallesCaldera_Mix.mat')
    case 3
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Hyytiala\Res_Hyytiala.mat')
    case 4
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Fort_Peck_Minnesota\Res_Fort_Peck.mat')
    case 5
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Morgan_Monroe\Res_MMSF.mat')
    case 6
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Scrub_Oak\Res_ScrubOak.mat')
    case 7
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Hainich\Res_Hainich.mat')
    case 8
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Chestnut_Ridge\Res_Chestnut_Ridge.mat')
    case 9
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Negrisia\Data_Negrisia.mat')
    case 10
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Ivotuk\Res_Ivotuk.mat')
    case 11
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Tharandt_site\Res_Tharandt.mat')
    case 12
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Niwot_Ridge_Forest\Res_Niwot_Rdige.mat')
    case 13
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\UKEbu\Res_EasterBush.mat')
    case 14
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Blodgett_Forest\Res_Blodgett.mat')
    case 15
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Castelporziano_site\Res_Castelporziano.mat')
    case 16
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Boreas_Black_Spruce\Res_BlackSpruce.mat')
    case 17
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\US_Ho1\Res_Howland_Main.mat')
    case 18
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\ESES1\Res_El_Saler.mat')
    case 19
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\FI_Sod\Res_Sodankyla.mat')
    case 20
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\BE_Vie\Res_Vielsalm.mat')
    case 21
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\FR_Hes\Res_Hesse.mat')
    case 22
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\DK_Sor\Res_Soroe.mat')
    case 23
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\US_Wcr\Res_Willow_Creek.mat')
    case 24
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\US_Bar\Res_Bartlett.mat')
    case 25
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Hyytiala\Res_Hyytiala.mat')
    case 26
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\FR_LBr\Res_LeBray.mat')
    case 27
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\NL_Loo\Res_Loobos.mat')
    case 28
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\RU_Fyo\Res_Fyodorovskoje.mat')
    case 29
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\US_Hem\Res_Hemlock.mat')
    case 30
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\US_LPH\Res_LPH.mat')
    case 31
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\FR_PUE\Res_Puechabon.mat')
    case 32
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\BE_Bra\Res_Braschaat.mat')
        Rn=0;
    case 33
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\ILYat\Res_Yatir.mat')
    case 34
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Silas\Res_Silas.mat')
    case 35
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\OPAL\Res_files\Res_BKS.mat')
    case 36
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\OPAL\Res_files\Res_PSO.mat')
    case 37
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\IEDri\Res_Dripsey.mat')
    case 38
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\ES_Amo\Res_Amoladeras.mat')
    case 39
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\ES_Agu\Res_Balsablanca.mat')
    case 40
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\ValldAlinya\Res_ValldAlinya.mat')
    case 41
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\IT_Noe\Res_Noe.mat')
    case 42
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Matra\Res_Matra.mat')
    case 43
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Ca_Obs\Res_SK_OBS.mat')
    case 44
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\ESLma\Res_EsLma.mat')
    case 45
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Metolius_Interm\Res_Metolius_Int.mat')
    case 46
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\ParkFalls\Res_Parkfalls.mat')
    case 47
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\ES_LJu\Res_EsLju.mat')
    case 48
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\ESES1\Res_El_Saler.mat')
    case 49
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Cedar_Bridge\Res_Cedar_Bridge.mat')
    case 50
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\New_Sites\Res_Howard_Springs.mat')
    case 51
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\New_Sites\Res_Adelaide_River.mat')
    case 52
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\New_Sites\Res_Mongo.mat')
    case 53
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\New_Sites\Res_Dry_River.mat')
    case 54
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\New_Sites\Res_Sturt_Plain.mat')
    case 55
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\New_Sites\Res_Sardinilla_Pasture.mat')
    case 56
        load('C:\DD\DESKTOP_SF\BILANCIO IDROLOGICO\TeCgam\Tibetan\InputData\Res_TibetanPE_New.mat')
    case 57
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Beringer_Flux\Res_Wallaby_New.mat');
        LE=LE2; H=H2; G=G2; Rn=Rn2;
    case 58
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Maui_Sugarcane_Windy\Res_Maui_Sugarcane_Windy.mat')
    case 59
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_RU-Ha1_FLUXNET2015_SUBSET_2002-2004_1-3\Res_Hakasia_Steppe.mat')
    case 60
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_AU-DaS_FLUXNET2015_FULLSET_2008-2014_2-3\Res_Daly_Savanna.mat')
    case 61
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_AU-DaP_FLUXNET2015_FULLSET_2007-2013_2-3\Res_Daly_Pasture.mat')
    case 62
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_AU-ASM_FLUXNET2015_SUBSET_2010-2014_2-3\Res_Alice_Springs.mat')
    case 63
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_AR-Vir_FLUXNET2015_FULLSET_2009-2012_1-3\Res_Virasoro.mat')
    case 64
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_AU-TTE_FLUXNET2015_FULLSET_2012-2013_1-3\Res_TiTree.mat')
    case 65
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_CN-Sw2_FLUXNET2015_FULLSET_2010-2012_1-3\Res_Siziwang_Grazed.mat')
    case 66
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_IT-Col_FLUXNET2015_SUBSET_1996-2014_1-3\Res_Collelongo.mat')
    case 67
        load('C:\DD\DESKTOP_SF\BILANCIO IDROLOGICO\Experimental Watershed\Caatinga\Res_Petrolina.mat')
        GPP1=GPP_f; GPP2=GPP_DT; GPP3=GPP1;GPP4=GPP1; H=H_f; LE=LE_f; Tsoil1=Ta*NaN;   SWC1 = Ta*0; Rn=Latm'-Lw_out+Rsw-Rsw_out;
    case 69
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_CH-Dav_FLUXNET2015_SUBSET_1997-2014_1-3\Res_Davos_New.mat')
    case 70
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Lavarone_Torgnon\Res_Lavarone_combined.mat')
        I=isnan(LE); GPP1(I)=NaN; GPP2(I)=NaN; GPP3(I)=NaN;  GPP4(I)=NaN;
    case 71
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Lavarone_Torgnon\Res_Torgnon.mat')
        GPP1(1:4180)=NaN; GPP2(1:4180)=NaN; GPP3(1:4180)=NaN;  GPP4(1:4180)=NaN;
    case 72
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_GF-Guy_FLUXNET2015_SUBSET_2004-2014_2-3\Res_Guyaflux.mat')
    case 73
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_IT-MBo_FLUXNET2015_SUBSET_2003-2013_1-3\Res_Monte_Bondone.mat')
    case 74
        load('C:\DD\DESKTOP_SF\BILANCIO IDROLOGICO\Experimental Watershed\Garmisch_Tereno\Res_Fendt.mat')
    case 75
        load('C:\DD\DESKTOP_SF\BILANCIO IDROLOGICO\Experimental Watershed\Garmisch_Tereno\Res_Rottenbuch.mat')
    case 76
        load('C:\DD\DESKTOP_SF\BILANCIO IDROLOGICO\Experimental Watershed\Garmisch_Tereno\Res_Graswang.mat')
    case 77
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_AU-Tum_FLUXNET2015_FULLSET_2001-2014_2-3\Res_Tumbarumba.mat')
    case 78
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_CN-Din_FLUXNET2015_SUBSET_2003-2005_1-3\Res_Dinghushan.mat')
    case 79
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_DK-ZaH_FLUXNET2015_FULLSET_2000-2014_2-3\Res_Zackenberg_Heath.mat')
    case 80
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_CN-Du2_FLUXNET2015_FULLSET_2006-2008_1-3\Res_Duolun_grassland_D01.mat')
        NEE_st=NEE1; NEE_or=NEE2;
    case 81
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_DE-Spw_FLUXNET2015_SUBSET_2010-2014_1-3\Res_Spreewald.mat')
    case 82
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Other_Sites\OPAL\Res_files\Res_PDF.mat')
    case 83
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_CN-Dan_FLUXNET2015_SUBSET_2004-2005_1-3\Res_Dangxiong.mat')
    case 84
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_RU-SkP_FLUXNET2015_SUBSET_2012-2014_1-4\Res_Yakutsk_Spasskaya.mat')
    case 85
        load('C:\DD\DESKTOP_SF\BILANCIO IDROLOGICO\Experimental Watershed\Greece_Tower\Res_NGreece.mat')
        GPP1=GPP; GPP2=GPP; GPP3=GPP;GPP4=GPP;  NEE_st=NEE; NEE_or=NEE;
    case 86
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_US-Ne1_FLUXNET2015_SUBSET_2001-2013_1-3\Res_US_Ne1.mat')
    case 87
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_US-Ne2_FLUXNET2015_SUBSET_2001-2013_1-3\Res_US_Ne2.mat')
    case 88
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_US-Ne3_FLUXNET2015_SUBSET_2001-2013_1-3\Res_US_Ne3.mat')
    case 89
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\AMF_US-Bo1_BASE-BADM_2-1\Res_Bondville1.mat')
        GPP1=GPP; GPP2=GPP; GPP3=GPP;GPP4=GPP;  NEE_st=NEE; NEE_or=NEE;
    case 90
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_US-Tw1_FLUXNET2015_SUBSET_2012-2014_1-4\Res_Twitchell_Wetland_West_Pond.mat')
        GPP1(1:4600)=NaN; GPP2(1:4600)=NaN; GPP3(1:4600)=NaN; GPP4(1:4600)=NaN;
    case 91
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_US-Twt_FLUXNET2015_SUBSET_2009-2014_1-3\Res_Twitchell_Island.mat')
    case 92
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_GL-ZaF_FLUXNET2015_SUBSET_2008-2011_2-4\Res_Zackenberg_Fen.mat')
    case 93
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_PH-RiF_FLUXNET-CH4_2012-2014_1-1\Res_PH_RiF.mat')
        GPP1(1:8135)=NaN; GPP2(1:8135)=NaN; GPP3(1:8135)=NaN; GPP4(1:8135)=NaN;
    case 94
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_MY-MLM_FLUXNET-CH4_2014-2015_1-1\Res_Maludam_National_Park.mat')
    case 95
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_CN-Qia_FLUXNET2015_SUBSET_2003-2005_1-4\Res_Qianyanzhou.mat')
    case 96
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_ID-Pag_FLUXNET-CH4_2016-2017_1-1\Res_Palangkaraya_undrained.mat')
    case 97
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_CN-Cng_FLUXNET2015_SUBSET_2007-2010_1-3\Res_Changling.mat')
    case 98
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_CN-Ha2_FLUXNET2015_SUBSET_2003-2005_1-3\Res_Haibei_Shrubland.mat')
    case 99
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_JP-SMF_FLUXNET2015_SUBSET_2002-2006_1-3\Res_Seto_Mixed_Forest.mat')
        Rn=Rsw;
    case 100
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_HK-MPM_FLUXNET-CH4_2016-2018_1-1\Res_Mai_Po_Mangrove.mat')
    case 101
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\AMF_US-MOz_BASE-BADM_11-5\Res_Missouri_Ozark.mat')
        GPP1=GPP; GPP2=GPP; GPP3=GPP;GPP4=GPP;  NEE_st=NEE; NEE_or=NEE;
    case 102
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\AMF_US-Skr_BASE-BADM_1-1\Res_US_Skr.mat')
        GPP1=GPP; GPP2=GPP; GPP3=GPP;GPP4=GPP;  NEE_st=NEE; NEE_or=NEE;
    case 103
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_JP-MBF_FLUXNET2015_SUBSET_2003-2005_1-3\Res_Moshiri_Birch_Forest.mat')
        Rn=Rsw;
    case 104
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_CN-HaM_FLUXNET2015_FULLSET_2002-2004_1-4\Res_HaibeiAlpineTibetsite.mat')
    case 105
        load('C:\DD\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\FLX_IT-Ro2_FLUXNET2015_SUBSET_2002-2012_1-3\Res_Roccarespampani2.mat')
        NEE_st=NEE2; NEE_or=NEE4;
end


%%%%%%%%%%%%%%%%
cd(cur_dir2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hs=H; Rns=Rn;
clear Lon SD SB Zbas DeltaGMT Lat
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
switch OPT_CASE
    case {1}
        load('Ris_Renon.mat')
        x1=1;
        x2=140256;
        yyy=1998:2013;
        Tsoil2 = Tsoil1;
        Idp=[ 50 150];
    case 2
        load('Ris_ValCalMix.mat')
        x1=1;
        x2=35064;
        yyy=2007:2010;
        Idp=[ 25 50];
    case 3
        load('Ris_HYT.mat');
        x1=1;
        x2=  140256;%
        yyy=1997:2012;
        Tsoil2=Ts2; Tsoil1=Ts1;
        Idp=[ 50 150];
    case 4
        load('Ris_FortPeck.mat');
        x1=1;
        x2= 61368;  %
        yyy=2000:2006;
        Idp=[ 100 200];
    case 5
        load('Ris_MMSF.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2= 70128;  %
        yyy=1999:2006;
    case 6
        load('Ris_Scrub_Oak.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=61368;  %
        yyy=2000:2006;
    case 7
        load('Ris_Hainich.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=70128;  %
        yyy=2000:2007;
    case 8
        load('Ris_Chestnut.mat');
        x1=1;
        x2=49464;  %
        yyy=2006:2010;
    case 9
        load('Ris_Negrisia.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=35064;  %
        yyy=2006:2009;
    case 10
        load('Ris_Ivotuk.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=28512;  %
        yyy=2004:2006;
    case 11
        load('Ris_Tharandt.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=122712;  %
        yyy=1997:2010;
    case 12
        load('Ris_NiwotRidge.mat');
        x1=1;
        x2=132960;  %
        yyy=1999:2013;
    case 13
        load('Ris_EasterBush.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=61368;  %
        yyy=2004:2010;
    case 14
        load('Ris_Blodgett.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=66788;  %
        yyy=2000:2006;
    case 15
        load('Ris_CastelPorz.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=78912;  %
        yyy=2000:2008;
    case 16
        load('Ris_BlackSpruce.mat');
        x1=1;
        x2=131490;  %
        yyy=1994:2008;
    case 17
        load('Ris_Howland1.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=70128;  %
        yyy=1997:2004;
    case 18
        load('Ris_ElSaler.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=70128;  %
        yyy=1999:2006;
    case 19
        load('Ris_Sodankyla.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=78912;  %
        yyy=2000:2008;
    case 20
        load('Ris_Vielsalm.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=87648;  %
        yyy=1997:2006;
    case 21
        load('Ris_Hesse.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=87648;  %
        yyy=1997:2006;
    case 22
        load('Ris_Soroe.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=87648;  %
        yyy=1997:2006;
    case 23
        load('Ris_WilCre.mat');
        %Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=140256;  %
        yyy=1998:2014;
    case 24
        load('Ris_Bartlett.mat');
        %Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=61368;  %
        yyy=2004:2010;
    case 25
        load('Ris_HYT.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=140256;  %
        yyy=1997:2012;
    case 26
        load('Ris_LeBray.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=70128;  %
        yyy=2000:2008;
    case 27
        load('Ris_Loobos.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=87648;  %
        yyy=1997:2006;
    case 28
        load('Ris_Fyodorovskoye.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=113952;  %
        yyy=1998:2010;
    case 29
        load('Ris_Ha2.mat');
        %Tsoil2=Ts2; Tsoil1=Ts1;
        Rns=Rn;
        x1=1;
        x2=52584;
        yyy=2005:2010;
    case 30
        load('Ris_LPH.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=70128;
        yyy=2003:2010;
    case 31
        load('Ris_Puechbon.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=105168;
        yyy=2000:2011;
    case 32
        load('Ris_Braaschat.mat');
        Tsoil2=Ts_2; Tsoil1=Ts_1;
        Rns=Rn; SWC1 = 0*Rns;
        x1=1;
        x2=157800;
        yyy=1996:2013;
    case 33
        load('Ris_Yatir.mat');
        Tsoil2=Ts2; Tsoil1=Ts1;
        x1=1;
        x2=78888;
        yyy=2001:2009;
    case 34
        load('Ris_Silas.mat');
        Tsoil2=Tsoil1; % Tsoil1=Ts1;
        x1=1;
        x2=87648;
        yyy=2005:2014;
    case 35
        load('Ris_BKS.mat');
        x1=1;
        x2=17520;
        yyy=2001:2002;
        GPP1 = Date*NaN; GPP2=GPP1; GPP3=GPP1; GPP4=GPP1; Tsoil2=Tsoil1;
    case 36
        load('Ris_PSO.mat');
        x1=1;
        x2=61368;
        yyy=2003:2009;
        GPP1 = Date*NaN; GPP2=GPP1; GPP3=GPP1; GPP4=GPP1; Tsoil2=Tsoil1;
    case 37
        load('Ris_Dripsey.mat');
        Tsoil1=Ts1; Tsoil2=Ts2;
        x1=1;
        x2=41672;
        yyy=2003:2007;
    case 38
        load('Ris_Amoladeras.mat')
        Tsoil2=Ts2;  Tsoil1=Ts1;
        x1=1;
        x2=43824;
        yyy=2007:2011;
    case 39
        load('Ris_Balsablanca.mat')
        Tsoil2=Ts2;  Tsoil1=Ts1;
        x1=1;
        x2=48960;
        yyy=2006:2011;
    case 40
        load('Ris_Valldalynia.mat')
        Tsoil2=Ts2;  Tsoil1=Ts1;
        x1=1;
        x2=43848;
        yyy=2004:2008;
    case 41
        load('Ris_Noe.mat')
        Tsoil2=Ts1;  Tsoil1=Ts1;
        x1=1;
        x2=50471;
        yyy=2004:2009;
    case 42
        load('Ris_Matra.mat')
        Tsoil2=Ts1;  Tsoil1=Ts1;
        x1=1;
        x2=43848;
        yyy=2004:2008;
    case 43
        load('Ris_SOBS.mat')
        Tsoil1 = SoilTemp_5cm;  Tsoil2=SoilTemp_10cm; GPP1 = GEP; GPP2=GEP; GPP3=GEP; GPP4 = GEP;
        NEE = NEP; % NEE2=NEE;
        SWC1 = NaN*GEP;
        x1=1;
        x2=166536;
        yyy=1997:2015;
    case 44
        load('Ris_EsLma.mat')
        Tsoil2=Ts1;  Tsoil1=Ts1;
        x1=1;
        x2=70128;
        yyy=2004:2011;
    case 45
        load('Ris_Metolius.mat')
        Tsoil2=Tsoil1;
        x1=1;
        x2=96432;
        yyy=2002:2012;
    case 46
        load('Ris_Parkfalls.mat')
        Rns=Rn;
        Tsoil2=Tsoil1;
        x1=1;
        x2=184080;
        yyy=1995:2015;
    case 47
        load('Ris_LJU.mat')
        Tsoil2=Ts1;  Tsoil1=Ts1;
        x1=1;
        x2=67247;
        yyy=2004:2011;
    case 48
        load('Ris_ElSaler.mat')
        Tsoil2=Ts1;  Tsoil1=Ts1;
        x1=1;
        x2=70128;
        yyy=1999:2006;
    case 49
        load('Ris_CedarBridge.mat')
        Tsoil2=Tsoil1;
        x1=1;
        x2=83304;
        yyy=1999:2006;
    case 50
        load('Ris_Howard_Springs.mat')
        Tsoil2=Tsoil1;  GPP1 = GPP; GPP3 = GPP; GPP4 =GPP;
        Hs(1:6000)=NaN; LE(1:6000)=NaN; Rns(1:6000)=NaN;
        x1=1;
        x2=122712;
        yyy=2001:2014;
        Idp=[ 150 250];
    case 51
        load('Ris_Adelaide_River.mat')
        Tsoil2=Tsoil1;
        LE(isnan(Rns))=NaN; Hs(isnan(Rns))=NaN;  GPP(isnan(Rns))=NaN; GPP2(isnan(Rns))=NaN;
        GPP1 = GPP; GPP3 = GPP; GPP4 =GPP;
        x1=1;
        x2=26304;
        yyy=2007:2009;
        Idp=[ 150 250];
    case 52
        load('Ris_Mongu.mat')
        Tsoil2=Tsoil1;  GPP1 = GPP; GPP3 = GPP; GPP4 =GPP;
        Hs(1:67610)=NaN; LE(1:67610)=NaN; Hs(83750:end)=NaN; LE(83750:end)=NaN;
        x1=1;
        x2=87672;
        yyy=2000:2009;
    case 53
        load('Ris_Dry_River.mat')
        Tsoil2=Tsoil1;  GPP1 = GPP; GPP3 = GPP; GPP4 =GPP;
        Hs(1:13000)=NaN; LE(1:13000)=NaN; Hs(47650:52600)=NaN; LE(47650:52600)=NaN;
        x1=1;
        x2=61368;
        yyy=2008:2014;
        Idp=[ 150 250];
    case 54
        load('Ris_Sturt_Plain.mat')
        Tsoil2=Tsoil1;  GPP1 = GPP; GPP3 = GPP; GPP4 =GPP;
        Hs(9800:11200)=NaN; LE(9800:11200)=NaN; Hs(1:5200)=NaN; LE(1:5200)=NaN;
        x1=1;
        x2=61368;
        yyy=2008:2014;
        Idp=[ 150 250];
    case 55
        load('Ris_Sardinilla.mat')
        Tsoil2=Tsoil1;  GPP1 = GPP; GPP3 = GPP; GPP4 =GPP;
        Hs(1:1950)=NaN; LE(1:1950)=NaN;
        x1=1;
        x2=26304;
        yyy=2007:2009;
        Idp=[ 150 250];
    case 56
        load('Ris_Tibetan_PE_bg.mat')
        GPP = NaN*Rn; Tsoil1 = NaN*Rn; Tsoil2=NaN*Rn;  GPP1 = GPP; GPP2 = GPP; GPP3 = GPP; GPP4 =GPP; SWC1 = NaN*Rn;
        x1=1;
        x2=20896;
        yyy=2015:2017;
    case 57
        load('Ris_Wallaby_FromSeed.mat')
        Tsoil1=Ts1; Tsoil2 = Ts1; SWC1=SWC1*100;
        GPP3 = GPP1; GPP4 =GPP1;
        x1=1;
        x2=162145;
        yyy=2000:2018;
        Idp=[ 150 250];
    case 58
        load('Ris_Maui_Sugarcane.mat')
        Ccrown=1.0;
        Tsoil2 = Tsoil1;
        GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=26304;
        yyy=2011:2013;
        Idp=[ 150 250];
    case 59
        load('Ris_Hakasia.mat')
        Tsoil2 = Tsoil1;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=26304;
        yyy=2002:2004;
        Idp=[ 150 250];
    case 60
        load('Ris_Daly_River_Savanna.mat')
        Tsoil2 = Tsoil1;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=61368;
        yyy=2008:2014;
        Idp=[ 150 250];
    case 61
        load('Ris_Daly_River_Pasture.mat')
        Tsoil2 = Tsoil1;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=61368;
        yyy=2007:2013;
        Idp=[ 150 250];
    case 62
        load('Ris_AliceSprings.mat')
        Tsoil2 = Tsoil1;  GPP1=GPP1*NaN; GPP2=GPP2*NaN;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=43824;
        yyy=2010:2014;
        Idp=[ 150 250];
    case 63
        load('Ris_Virasoro.mat')
        Tsoil2 = Tsoil1;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=35064;
        yyy=2009:2012;
        Idp=[ 150 250];
    case 64
        load('Ris_TiTree.mat')
        Tsoil2 = Tsoil1;  GPP1=GPP1*NaN; GPP2=GPP2*NaN;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=17543;
        yyy=2012:2013;
        Idp=[ 150 250];
    case 65
        load('Ris_SiziwangGrazed.mat')
        Tsoil2 = Tsoil1; GPP1=GPP1*NaN; GPP2=GPP2*NaN;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=26304;
        yyy=2010:2012;
        Idp=[ 150 250];
    case 66
        load('Ris_Collelongo.mat')
        Tsoil2 = Tsoil1;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=166560;
        yyy=1996:2014;
        Idp=[ 150 250];
    case 67
        load('Ris_Petrolina_New2.mat')
        Tsoil2 = Tsoil1;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=70127;
        yyy=2011:2018;
        Idp=[ 150 250];
    case 69
        load('Ris_Davos.mat')
        Tsoil2 = Tsoil1;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=157776;
        yyy=1997:2014;
        Idp=[ 150 250];
    case 70
        load('Ris_Lavarone.mat')
        Tsoil2 = Tsoil1;
        x1=1;
        x2=131496;
        yyy=2000:2014;
        Idp=[ 150 250];
    case 71
        load('Ris_Torgnon.mat')
        Tsoil2 = Tsoil1;
        x1=1;
        x2=61368;
        yyy=2008:2014;
        Idp=[ 150 250];
    case 72
        load('Ris_Guyaflux.mat')
        Tsoil2 = Tsoil1;
        x1=1;
        x2=96432;
        yyy=2004:2014;
        Idp=[ 150 250];
    case 73
        load('Ris_MonteBondone.mat')
        Tsoil2 = Tsoil1;
        x1=1;
        x2=96432;
        yyy=2003:2013;
        Idp=[ 150 250];
    case 74
        load('Ris_Fendt.mat')
        %Depth 50 35 25 12 6 2 cm
        GPP3 = GPP1; GPP4 =GPP1; GPP2 =GPP1; GPP1 = GPP1;
        Tsoil2 =Tsoil(:,2);
        Tsoil1 =Tsoil(:,4);
        SWC1 = SWC(:,4);
        SWC2 = SWC(:,2);
        x1=1;
        x2=61368;
        yyy=2012:2018;
        Idp=[ 120 350];
    case 75
        load('Ris_Rottenbuch.mat')
        %Depth 50 35 25 12 6 2 cm
        GPP3 = GPP1; GPP4 =GPP1; GPP2 =GPP1; GPP1 = GPP1;
        Tsoil2 =Tsoil(:,2);
        Tsoil1 =Tsoil(:,4);
        SWC1 = SWC(:,4);
        SWC2 = SWC(:,2);
        x1=1;
        x2=61368;
        yyy=2012:2018;
        Idp=[ 120 350];
    case 76
        load('Ris_Graswang.mat')
        %Depth 50 35 25 12 6 2 cm
        GPP3 = GPP1; GPP4 =GPP1; GPP2 =GPP1; GPP1 = GPP1;
        Tsoil2 =Tsoil(:,2);
        Tsoil1 =Tsoil(:,4);
        SWC1 = SWC(:,4);
        SWC2 = SWC(:,2);
        x1=1;
        x2=61368;
        yyy=2012:2018;
        Idp=[ 120 350];
    case 77
        load('Ris_Tumbarumba.mat')
        Tsoil2 = Tsoil1;
        GPP1=GPP3; GPP2=GPP4;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=122712;
        yyy=2001:2014;
        Idp=[ 150 250];
    case 78
        load('Ris_Dinghushan.mat')
        Tsoil2 = Tsoil1;
        %GPP1=GPP3; GPP2=GPP4;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=26304;
        yyy=2003:2005;
        Idp=[ 150 250];
    case 79
        load('Ris_Zackenberg_Heath.mat')
        Tsoil2 = Tsoil1;
        %GPP1=GPP3; GPP2=GPP4;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=131496;
        yyy=2000:2014;
        Idp=[50 100];
    case 80
        load('Ris_Duolun.mat')
        Tsoil2 = Tsoil1; % GPP1=GPP1*NaN; GPP2=GPP2*NaN;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=26304;
        yyy=2006:2008;
        Idp=[ 150 250];
    case 81
        load('Ris_Spreewald.mat')
        Tsoil2 = Tsoil1; %
        x1=1;
        x2=43824;
        yyy=2010:2014;
        Idp=[ 150 250];
    case 82
        load('Ris_PDF.mat');
        x1=1;
        x2=35064;
        yyy=2002:2005;
        GPP1 = Date'*NaN; GPP2=GPP1; GPP3=GPP1; GPP4=GPP1; Tsoil2=Tsoil1;
        SWC1=SWC1*100;
        Idp=[ 150 250];
    case 83
        load('Ris_Dangxiong.mat')
        Tsoil2 = Tsoil1; %
        x1=1;
        x2=17544;
        yyy=2004:2005;
        Idp=[ 150 250];
    case 84
        load('Ris_Yakutsk.mat')
        Tsoil2 = Tsoil1; %
        x1=1;
        x2=26304;
        yyy=2012:2014;
        Idp=[ 150 250];
    case 85
        load('Ris_NGreece.mat')
        Tsoil1 = Tsoil; %
        Tsoil2 = Tsoil; %
        SWC1=SWC10cm;
        x1=1;
        x2=14568;
        yyy=2019:2021;
        Idp=[ 100 150];
    case 86
        load('Ris_Mead1.mat')
        Ccrown=[1.0];
        x1=1;
        x2=113952;
        yyy=2001:2013;
        Idp=[ 150 250];
    case 87
        load('Ris_Mead2.mat')
        Ccrown=[1.0 1.0];
        x1=1;
        x2=113952;
        yyy=2001:2013;
        Idp=[ 150 250];
    case 88
        load('Ris_Mead3.mat')
        Ccrown=[1.0 1.0];
        x1=1;
        x2=113952;
        yyy=2001:2013;
        Idp=[ 150 250];
    case 89
        load('Ris_Bondville1.mat')
        Ccrown=[1.0 1.0];
        x1=1;
        x2=113976;
        yyy=1996:2008;
        Idp=[ 150 250];
    case 90
        load('Ris_Tw_WestPond.mat')
        x1=1;
        x2=26304;
        yyy=2012:2014;
        Idp=[ 150 250];
    case 91
        load('Ris_Tw_Rice.mat')
        Ccrown=1.0;
        x1=1;
        x2=52584;
        yyy=2009:2014;
        Idp=[ 150 250];
    case 92
        load('Ris_Zackenberg_Fen.mat')
        Tsoil2 = Tsoil1;
        %GPP1=GPP3; GPP2=GPP4;
        %GPP3 = GPP; GPP4 =GPP; GPP2 =GPP; GPP1 = GPP;
        x1=1;
        x2=35064;
        yyy=2008:2011;
        Idp=[50 100];
    case 93
        load('Ris_PH_Rif.mat')
        Ccrown=1.0;
        x1=1;
        x2=26304;
        yyy=2012:2014;
        Idp=[ 150 250];
    case 94
        load('Ris_Maludam.mat')
        Rns=Rn;
        x1=1;
        x2=17519;
        yyy=2014:2015;
        Idp=[ 150 250];
    case 95
        load('Ris_Qianyanzhou.mat')
        x1=1;
        x2=26304;
        yyy=2003:2005;
        Idp=[ 150 250];
    case 96
        load('Ris_ID_Pag.mat')
        x1=1;
        x2=8759;
        yyy=2016:2017;
        Idp=[ 150 250];
    case 97
        load('Ris_Changling.mat')
        x1=1;
        x2=35064;
        yyy=2007:2010;
        Idp=[ 150 250];
    case 98
        load('Ris_Haibei_2.mat')
        x1=1;
        x2=26304;
        yyy=2003:2005;
        Idp=[ 150 250];
    case 99
        load('Ris_Seto_Mixed_Forest.mat')
        x1=1;
        x2=43824;
        yyy=2002:2006;
        Idp=[ 150 250];
    case 100
        load('Ris_Mai_Po.mat')
        x1=1;
        x2=26303;
        yyy=2016:2018;
        Idp=[ 150 250];
    case 101
        load('Ris_Missouri_Ozark.mat')
        x1=1;
        x2=148775;
        yyy=2005:2021;
        Idp=[ 150 250];
    case 102
        load('Ris_US_Skr.mat')
        x1=1;
        x2=43824;
        yyy=2007:2011;
        Idp=[ 150 250];
    case 103
        load('Ris_Moshiri_Birch_Forest.mat')
        x1=1;
        x2=26304;
        yyy=2003:2005;
        Idp=[ 150 250];
    case 104
        load('Ris_Haibei_Meadow.mat')
        x1=1;
        x2=26304;
        yyy=2002:2004;
        Idp=[ 150 250];
    case 105
        load('Ris_Ro2_bg.mat')
        x1=1;
        x2=96432;
        yyy=2002:2012;
        Idp=[ 150 250];
end
yym1 = yyy(1)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zs= cumsum(Dz);
Rns=Rns(x1:x2);
LE=LE(x1:x2);
Hs=Hs(x1:x2);
GPP1=GPP1(x1:x2);
GPP2=GPP2(x1:x2);
GPP3=GPP3(x1:x2);
GPP4=GPP4(x1:x2);
SWC1=SWC1(x1:x2);
Ts2=Tsoil2(x1:x2);
Ts1=Tsoil1(x1:x2);
%%%%%%%%%%%%%%%%
RnV =dQVEG +HV+QEV+Lpho.*(TsV~=0);
H=H(x1:x2);
Rn=Rn(x1:x2);
QE=QE(x1:x2);
Datam=Datam(x1:x2,:);
An_L=An_L(x1:x2,:);
Rdark_L =Rdark_L(x1:x2,:);
An_H=An_H(x1:x2,:);
Rdark_H =Rdark_H(x1:x2,:);
%Date=Date(x1:x2);
%%%%%%%%%%%%%%%%
OPT_SoilBiogeochemistry=1;
%%%%
OPT_COMP_VAL=1;
%%%
OPT_SN=0;
%%%%
if OPT_SN
    Hs(Csno==1)=NaN;
    Rns(Csno==1)=NaN;
    LE(Csno==1)=NaN;
    %%%%%
    Hs(TsV~=0)=NaN;
    Rns(TsV~=0)=NaN;
    LE(TsV~=0)=NaN;
else
    Rn=Rn+RnV;
    H=H+HV;
    QE=QE+QEV;
end
%%%%%%%%%%%%%%%%% GPP analysis
GPPs = [ GPP1 GPP2 GPP3 GPP4];
GPPs =nanmean(GPPs,2);
GPPs(GPPs<0)=0;
%%%%%%%%%%%%%%%%%%%%%%
GPP =(An_H+Rdark_H+An_L+Rdark_L)*Ccrown';
fr=24;
m=floor(length(GPPs)/fr);

Xm=reshape(GPPs(1:m*fr),fr,m);
Xs=reshape(GPP(1:m*fr),fr,m);
GPPd_m=1.0368*mean(Xm);
GPPd_s=1.0368*mean(Xs);

if OPT_CASE ==46; GPPd_m=1.0368*nanmean(Xm); end

Dd=Date(1):1:Date(end);
if length(Dd) > length(GPPd_m)
    Dd=Dd(1:end-1);
end
jDay=jDay(1:length(Dd));

if OPT_SoilBiogeochemistry == 1
    NEEs = 0.5*(NEE_st+NEE_or);
    Reco2 = 0.5*(Reco2+Reco1);
    %%%%%%%%%%%%%%%%%
    NEE = -(NPP_H+NPP_L)*Ccrown' + R_litter + R_microbe +R_ew; %% [gC/m2 day]
    RE=  R_litter + R_microbe +R_ew +(RA_H+RA_L)*Ccrown' ;%% [gC/m2 day]
    NEE=NEE(2:end); RE=RE(2:end);
    fr=24;
    m=floor(length(NEEs)/fr);
    Xm=reshape(NEEs(1:m*fr),fr,m);
    NEEd=1.0368*mean(Xm);%% [gC/m2 day]
    NEEd_s = NEE;
    %%%%
    Xm=reshape(Reco2(1:m*fr),fr,m);
    REd=1.0368*mean(Xm);%% [gC/m2 day]
    %%%%
    REd_s = RE;
    %%%%
end
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% FLUXES
figure(1)
set(gca,'fontsize',12,'fontweight','normal');
subplot(3,1,2)
plot(Hs)
hold on ; grid on ;
plot(H,'--r')
title('Sensible Heat')
legend('OBS','SIM')
xlabel('Hours'); ylabel('[W m^{-2}')
subplot(3,1,1)
plot(Rns)
hold on ; grid on ;
plot(Rn,'--r')
title('Net Radiation')
xlabel('Hours'); ylabel('[W m^{-2}')
legend('OBS','SIM')
subplot(3,1,3)
plot(LE)
hold on ; grid on ;
plot(QE,'--r')
title('Latent Heat')
xlabel('Hours'); ylabel('[W m^{-2}')
legend('OBS','SIM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SCATTER PLOTS
figure(2)
set(gca,'fontsize',12,'fontweight','normal');
set(gcf,'color','w');
%
de=10;
edgs=[-200:de:1000];
edges{1}=edgs-de/2;
edges{2}=edgs+de/2;
colormap(jet); %parula; $jet
cmap = colormap;
cmap(1,:) = [0.98 0.98 0.98];
colormap(cmap);
mnp = 100;

subplot(1,3,1)
[values, centers] = hist3([Rn Rns],'Edges',edges);
imagesc(centers{:}, values')
axis xy
colorbar
caxis([0 mnp])
%hold on
%plot(RIS_SITE.Rn,SITE.Rn,'kx')
xlabel('SIM R_n')
ylabel('OBS R_n')
R=corrcoef(Rn(not(isnan(Rns))),Rns(not(isnan(Rns))));
hold on ;
grid on ;
plot(edgs,edgs,'r')
title(strcat('a R =', num2str(R(1,2))),'FontSize',12,'FontWeight','normal')
xlim([-200 1000]) ; ylim([-200 1000])
%
subplot(1,3,2)
[values, centers] = hist3([H Hs],'Edges',edges);
imagesc(centers{:}, values')
axis xy
colorbar
caxis([0 mnp])
%hold on
%plot(RIS_SITE.H,SITE.H,'kx')
xlabel('SIM H')
ylabel('OBS H')
R=corrcoef(H(not(isnan(Hs))),Hs(not(isnan(Hs))));
hold on ;
grid on ;
plot(edgs,edgs,'r')
title(strcat('b) R =', num2str(R(1,2))),'FontSize',12,'FontWeight','normal')
xlim([-200 800]) ; ylim([-200 800])
%
subplot(1,3,3)
[values, centers] = hist3([QE LE],'Edges',edges);
imagesc(centers{:}, values')
axis xy
colorbar
caxis([0 mnp])
%hold on
%plot(RIS_SITE.QE,SITE.LE,'kx')
xlabel('SIM \lambdaE')
ylabel('OBS \lambdaE')
R=corrcoef(QE(not(isnan(LE))),LE(not(isnan(LE))));
hold on ;
grid on ;
plot(edgs,edgs,'r')
title(strcat('c) R =', num2str(R(1,2))),'FontSize',12,'FontWeight','normal')
xlim([-200 600]) ; ylim([-200 600])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE ONLY WHERE THERE ARE VALUES
if OPT_COMP_VAL == 1
    Rn(isnan(Rns))=NaN;
    QE(isnan(LE))=NaN;
    H(isnan(Hs))=NaN;
end

%%%%% ANNUAL Rn
YY=year(Date);
for i=yyy
    I=year(Date)==i ;
    Rnym(i-yym1)=nanmean(Rns(I));
    Rnys(i-yym1)=nanmean(Rn(I));
end

%%%%% ANNUAL LE
YY=year(Date);
for i=yyy
    I=year(Date)==i ;
    LEym(i-yym1)=nanmean(LE(I));
    LEys(i-yym1)=nanmean(QE(I));
end

%%%%% ANNUAL H
YY=year(Date);
for i=yyy
    I=year(Date)==i ;
    Hym(i-yym1)=nanmean(Hs(I));
    Hys(i-yym1)=nanmean(H(I));
end

%%%%% ANNUAL GPP
YY=year(Dd);
for i=yyy
    I=year(Dd)==i;
    GPPym(i-yym1)=mean(GPPd_m(I));
    GPPys(i-yym1)=mean(GPPd_s(I));
end


figure(4)
subplot(3,1,1)
plot(yyy,Rnym(:),'-ob')
a = polyfit(yyy,Rnym(:)',1);
slm=a(1);
hold on
plot(yyy,Rnys(:),'-xr')
a = polyfit(yyy,Rnys(:)',1);
sls=a(1);
xlim([min(yyy)-1 max(yyy)+1]);
xlabel('Year'); ylabel('Rn [W m^{-2}]')
title(strcat('OBS Trend=',num2str(slm),'; SIM Trend=',num2str(sls)))
legend('OBS','SIM')
ylim([10 180])
subplot(3,1,2)
plot(yyy,LEym(:),'-ob')
a = polyfit(yyy,LEym(:)',1);
slm=a(1);
hold on
plot(yyy,LEys(:),'-xr')
a = polyfit(yyy,LEys(:)',1);
sls=a(1);
xlim([min(yyy)-1 max(yyy)+1]);
xlabel('Year'); ylabel('\lambda E [W m^{-2}]')
title(strcat('OBS Trend=',num2str(slm),'; SIM Trend=',num2str(sls)))
legend('OBS','SIM')
ylim([0 120])
subplot(3,1,3)
plot(yyy,Hym(:),'-ob')
a = polyfit(yyy,Hym(:)',1);
slm=a(1);
hold on
plot(yyy,Hys(:),'-xr')
a = polyfit(yyy,Hys(:)',1);
sls=a(1);
xlim([min(yyy)-1 max(yyy)+1]);
xlabel('Year'); ylabel('H [W m^{-2}]')
title(strcat('OBS Trend=',num2str(slm),'; SIM Trend=',num2str(sls)))
legend('OBS','SIM')
ylim([-20 120])


figure(114)
subplot(1,1,1)
plot(yyy,GPPym(:),'-ob')
a = polyfit(yyy,GPPym(:)',1);
slm=a(1);
hold on
plot(yyy,GPPys(:),'-xr')
a = polyfit(yyy,GPPys(:)',1);
sls=a(1);
xlim([min(yyy)-1 max(yyy)+1]);
xlabel('Year'); ylabel('[gC m^{-2} day^{-1}]')
title(strcat('OBS Trend=',num2str(slm),'; SIM Trend=',num2str(sls)))
legend('OBS','SIM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  SEASONAL FLUXES

%%%%% Daily Rn

Xm=reshape(Rns(1:m*fr),fr,m);
Xs=reshape(Rn(1:m*fr),fr,m);
Rnd_m=mean(Xm);
Rnd_s=mean(Xs);

for j=1:365
    for i=yyy
        I=year(Dd)==i &  jDay'==j ;
        Rn_doy(i-yym1)=mean(Rnd_m(I));
        Rn_doy2(i-yym1)=mean(Rnd_s(I));
    end
    M_Rn(j) = nanmean(Rn_doy);
    M_Rn_sim(j) = nanmean(Rn_doy2);
    clear  Rn_doy Rn_doy2
end

SLOPEp = M_Rn;
window_S =30;
SLO_F = 0*SLOPEp;
V=[ SLOPEp , SLOPEp , SLOPEp];
Vf= filter(ones(1,window_S )/window_S,1,V);
Vf(1:window_S )=Vf(window_S +1);
SLO_F=Vf([366:730]+window_S/2);
M_Rn_smooth = SLO_F;
%%%
SLOPEp = M_Rn_sim;
window_S =30;
SLO_F = 0*SLOPEp;
V=[ SLOPEp , SLOPEp , SLOPEp];
Vf= filter(ones(1,window_S )/window_S,1,V);
Vf(1:window_S )=Vf(window_S +1);
SLO_F=Vf([366:730]+window_S/2);
M_Rn_sim_smooth = SLO_F;

figure(5)
subplot(1,3,1)
plot(1:365,M_Rn,'ob')
hold on
plot(1:365,M_Rn_sim,'xr')
plot(1:365,M_Rn_smooth,'-b','Linewidth',2);
plot(1:365,M_Rn_sim_smooth,'-r','Linewidth',2);
xlabel('Doy'); ylabel('[W m^{-2}]')
legend('OBS','SIM')
title('Rn')

%%%%% Daily LE

Xm=reshape(LE(1:m*fr),fr,m);
Xs=reshape(QE(1:m*fr),fr,m);
LEd_m=mean(Xm);
LEd_s=mean(Xs);

%%% Trend for LE - Seasonal
for j=1:365
    for i=yyy
        I=year(Dd)==i &  jDay'==j ;
        LE_doy(i-yym1)=mean(LEd_m(I));
        LE_doy2(i-yym1)=mean(LEd_s(I));
    end
    M_LE(j) = nanmean(LE_doy);
    M_LE_sim(j) = nanmean(LE_doy2);
    clear  LE_doy LE_doy2
end

SLOPEp = M_LE;
window_S =30;
SLO_F = 0*SLOPEp;
V=[ SLOPEp , SLOPEp , SLOPEp];
Vf= filter(ones(1,window_S )/window_S,1,V);
Vf(1:window_S )=Vf(window_S +1);
SLO_F=Vf([366:730]+window_S/2);
M_LE_smooth = SLO_F;
%%%
SLOPEp = M_LE_sim;
window_S =30;
SLO_F = 0*SLOPEp;
V=[ SLOPEp , SLOPEp , SLOPEp];
Vf= filter(ones(1,window_S )/window_S,1,V);
Vf(1:window_S )=Vf(window_S +1);
SLO_F=Vf([366:730]+window_S/2);
M_LE_sim_smooth = SLO_F;

figure(5)
subplot(1,3,2)
plot(1:365,M_LE,'ob')
hold on
plot(1:365,M_LE_sim,'xr')
plot(1:365,M_LE_smooth,'-b','Linewidth',2);
plot(1:365,M_LE_sim_smooth,'-r','Linewidth',2);
xlabel('Doy'); ylabel('[W m^{-2}]')
legend('OBS','SIM')
title('\lambda E')

%%%%% Daily H

Xm=reshape(Hs(1:m*fr),fr,m);
Xs=reshape(H(1:m*fr),fr,m);
Hd_m=mean(Xm);
Hd_s=mean(Xs);

%%% Trend for H - Seasonal
for j=1:365
    for i=yyy
        I=year(Dd)==i &  jDay'==j ;
        H_doy(i-yym1)=mean(Hd_m(I));
        H_doy2(i-yym1)=mean(Hd_s(I));
    end
    M_H(j) = nanmean(H_doy);
    M_H_sim(j) = nanmean(H_doy2);
    clear  H_doy H_doy2
end

SLOPEp = M_H;
window_S =30;
SLO_F = 0*SLOPEp;
V=[ SLOPEp , SLOPEp , SLOPEp];
Vf= filter(ones(1,window_S )/window_S,1,V);
Vf(1:window_S )=Vf(window_S +1);
SLO_F=Vf([366:730]+window_S/2);
M_H_smooth = SLO_F;
%%%
SLOPEp = M_H_sim;
window_S =30;
SLO_F = 0*SLOPEp;
V=[ SLOPEp , SLOPEp , SLOPEp];
Vf= filter(ones(1,window_S )/window_S,1,V);
Vf(1:window_S )=Vf(window_S +1);
SLO_F=Vf([366:730]+window_S/2);
M_H_sim_smooth = SLO_F;

figure(5)
subplot(1,3,3)
plot(1:365,M_H,'ob')
hold on
plot(1:365,M_H_sim,'xr')
plot(1:365,M_H_smooth,'-b','Linewidth',2);
plot(1:365,M_H_sim_smooth,'-r','Linewidth',2);
xlabel('Doy'); ylabel('[W m^{-2}]')
legend('OBS','SIM')
title('H')



%%% Trend for GPP - Seasonal
for j=1:365
    for i=yyy
        I=year(Dd)==i &  jDay'==j ;
        GPP_doy(i-yym1)=mean(GPPd_m(I));
        GPP_doy2(i-yym1)=mean(GPPd_s(I));
    end
    M_GPP(j) = nanmean(GPP_doy);
    M_GPP_sim(j) = nanmean(GPP_doy2);
    clear  GPP_doy GPP_doy2
end

SLOPEp = M_GPP;
window_S =30;
SLO_F = 0*SLOPEp;
V=[ SLOPEp , SLOPEp , SLOPEp];
Vf= filter(ones(1,window_S )/window_S,1,V);
Vf(1:window_S )=Vf(window_S +1);
SLO_F=Vf([366:730]+window_S/2);
M_GPP_smooth = SLO_F;
%%%
SLOPEp = M_GPP_sim;
window_S =30;
SLO_F = 0*SLOPEp;
V=[ SLOPEp , SLOPEp , SLOPEp];
Vf= filter(ones(1,window_S )/window_S,1,V);
Vf(1:window_S )=Vf(window_S +1);
SLO_F=Vf([366:730]+window_S/2);
M_GPP_sim_smooth = SLO_F;


figure(115)
subplot(1,1,1)
plot(1:365,M_GPP,'ob')
hold on
plot(1:365,M_GPP_sim,'xr')
plot(1:365,M_GPP_smooth,'-b','Linewidth',2);
plot(1:365,M_GPP_sim_smooth,'-r','Linewidth',2);
xlabel('Doy'); ylabel('[gC m^{-2} day^{-1} year^{-1}]')
legend('OBS','SIM')
title('GPP')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% DAILY CYCLES

dt=1;
n=length(Rn); fr=24/dt;
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);
Xp=reshape(Rns(1:m*fr),fr,m);
Xsp=reshape(Rn(1:m*fr),fr,m);
Xm=nanmean(Xp');
Xs=nanmean(Xsp');
stdXm=nanstd(Xp');
stdXs=nanstd(Xsp');
%t=1:fr;
t=[Datam(1,4):23, 0:Datam(1,4)-1];
figure(6)
subplot(1,3,1)
set(gca,'Fontsize',11)
plot(t,Xm,'r','LineWidth', 2);
hold on ; grid on;
plot(t,Xs,'--b','LineWidth', 2);
plot(t,stdXm,'^r','LineWidth', 2);
plot(t,stdXs,'^b','LineWidth', 2);
title('NET RADIATION')
legend('OBS','SIM')
xlabel('Hour')
ylabel('[W m^{-2}]')
mR=max(max(Xm,Xs))+50;
axis([0 24 -100 mR])
clear Xm Xs Xp Xsp stdXm stdXs t mR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


m=floor(n/fr);
Xp=reshape(LE(1:m*fr),fr,m);
Xsp=reshape(QE(1:m*fr),fr,m);
Xm=nanmean(Xp');
Xs=nanmean(Xsp');
stdXm=nanstd(Xp');
stdXs=nanstd(Xsp');
%t=1:fr;
t=[Datam(1,4):23, 0:Datam(1,4)-1];
figure(6)
subplot(1,3,2);
set(gca,'Fontsize',11)
plot(t,Xm,'r','LineWidth', 2);
hold on ; grid on;
plot(t,Xs,'--b','LineWidth', 2);
plot(t,stdXm,'^r','LineWidth', 2);
plot(t,stdXs,'^b','LineWidth', 2);
title('LATENT HEAT')
legend('OBS','SIM')
xlabel('Hour')
ylabel('[W m^{-2}]')
mR=max(max(Xm,Xs))+50;
axis([0 24 -20 mR-30])
clear Xm Xs Xp Xsp stdXm stdXs t mR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


m=floor(n/fr);
Xp=reshape(Hs(1:m*fr),fr,m);
Xsp=reshape(H(1:m*fr),fr,m);
Xm=nanmean(Xp');
Xs=nanmean(Xsp');
stdXm=nanstd(Xp');
stdXs=nanstd(Xsp');
%t=1:fr;
t=[Datam(1,4):23, 0:Datam(1,4)-1];
figure(6)
subplot(1,3,3);
set(gca,'Fontsize',11)
plot(t,Xm,'r','LineWidth', 2);
hold on ; grid on;
plot(t,Xs,'--b','LineWidth', 2);
plot(t,stdXm,'^r','LineWidth', 2);
plot(t,stdXs,'^b','LineWidth', 2);
title('SENSIBLE HEAT')
legend('OBS','SIM')
xlabel('Hour')
ylabel('[W m^{-2}]')
mR=max(max(Xm,Xs))+50;
axis([0 24 -60 mR-30])

clear Xm Xs Xp Xsp stdXm stdXs t mR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Smoothed series


figure(7)

window_S =24*15;

V=Rns;  Vf=NaN*V;
for k=1:(length(V)-window_S)
    Vf(k) = nanmean(V(k:k+window_S));
end
Rns_smooth = Vf;
%%%
V=Rn; Vf=NaN*V;
for k=1:(length(V)-window_S)
    Vf(k) = nanmean(V(k:k+window_S));
end
Rn_sim_smooth = Vf;

subplot(2,2,1)
%h1=plot(Date,Rns,'--b','Linewidth',0.1);
%hold on
%h2=plot(Date,Rn,'--r','Linewidth',0.1);
h3=plot(Date,Rns_smooth,'-b','Linewidth',2);
hold on
h4=plot(Date,Rn_sim_smooth,'-r','Linewidth',2);
datetick('x',11)
xlabel('Date'); ylabel('[W m^{-2}]')
%legend([h3 h4],{'OBS','SIM'})
legend('OBS','SIM')
title('Rn')
grid on


V=Hs;  Vf=NaN*V;
for k=1:(length(V)-window_S)
    Vf(k) = nanmean(V(k:k+window_S));
end
Hs_smooth = Vf;
%%%
V=H; Vf=NaN*V;
for k=1:(length(V)-window_S)
    Vf(k) = nanmean(V(k:k+window_S));
end
H_sim_smooth = Vf;

subplot(2,2,2)
%h1=plot(Date,Hs,'--b','Linewidth',0.1);
%hold on
%h2=plot(Date,H,'--r','Linewidth',0.1);
h3=plot(Date,Hs_smooth,'-b','Linewidth',2);
hold on
h4=plot(Date,H_sim_smooth,'-r','Linewidth',2);
datetick('x',11)
xlabel('Date'); ylabel('[W m^{-2}]')
%legend([h3 h4],{'OBS','SIM'})
legend('OBS','SIM')
title('H')
grid on


V=LE;  Vf=NaN*V;
for k=1:(length(V)-window_S)
    Vf(k) = nanmean(V(k:k+window_S));
end
LE_smooth = Vf;
%%%
V=QE; Vf=NaN*V;
for k=1:(length(V)-window_S)
    Vf(k) = nanmean(V(k:k+window_S));
end
LE_sim_smooth = Vf;

subplot(2,2,3)
%h1=plot(Date,LE,'--b','Linewidth',0.1);
%hold on
%h2=plot(Date,QE,'--r','Linewidth',0.1);
h3=plot(Date,LE_smooth,'-b','Linewidth',2);
hold on
h4=plot(Date,LE_sim_smooth,'-r','Linewidth',2);
datetick('x',11)
xlabel('Date'); ylabel('[W m^{-2}]')
%legend([h3 h4],{'OBS','SIM'})
legend('OBS','SIM')
title('\lambdaE')
grid on

GPP =(An_H+Rdark_H+An_L+Rdark_L)*Ccrown';

V=GPPs;  Vf=NaN*V;
for k=1:(length(V)-window_S)
    Vf(k) = nanmean(V(k:k+window_S));
end
GPPs_smooth = Vf;
%%%
V=GPP; Vf=NaN*V;
for k=1:(length(V)-window_S)
    Vf(k) = nanmean(V(k:k+window_S));
end
GPP_sim_smooth = Vf;

subplot(2,2,4)
%h1=plot(Date,GPPs,'--b','Linewidth',0.1);
%hold on
%h2=plot(Date,GPP,'--r','Linewidth',0.1);
h3=plot(Date,GPPs_smooth,'-b','Linewidth',2);
hold on
h4=plot(Date,GPP_sim_smooth,'-r','Linewidth',2);
datetick('x',11)
xlabel('Date'); ylabel('[umol m^{-2} s^{-1}]')
%legend([h3 h4],{'OBS','SIM'})
legend('OBS','SIM')
title('GPP')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Old plots

figure(8)
n=length(Rn);
m=floor(n/fr);
Xp=reshape(An_H(1:m*fr)+An_L(1:m*fr),fr,m);
Xsp=reshape(Rdark_H(1:m*fr)+Rdark_L(1:m*fr),fr,m);
Xm=nanmean(Xp');
Xs=mean(Xsp');
%t=1:fr;
t=[Datam(1,4):23, 0:Datam(1,4)-1];
subplot(1,1,1);
set(gca,'FontSize',10);
plot(t,Xm,'b','LineWidth', 3);
hold on ; grid on;
plot(t,Xs,'--r','LineWidth', 3);
legend('Net Assimiliation Rate','Foliage Respiration')
ylabel('[\mu mol CO_2 / m^2 s ]'); xlabel('Hour');
clear Xm Xs Xp Xsp stdXm stdXs t mR

figure(9)
plot(cumsum(QE(not(isnan(LE)))),'b','LineWidth', 2);
hold on ; grid on
plot(cumsum(LE(not(isnan(LE)))),'--r','LineWidth', 2);
legend('SIM','OBS')
ylabel('Cumulative \lambda E'); xlabel('Date');


dt=1;
n=length(Rn); fr=24/dt;
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);
Xp=reshape(Rns(1:m*fr),fr,m);
Xsp=reshape(Rn(1:m*fr),fr,m);
Xm=nanmean(Xp');
Xs=nanmean(Xsp');
t=[Datam(1,4):23, 0:Datam(1,4)-1];
figure(10)
subplot(1,1,1);
set(gca,'FontSize',10);
plot(t,Xm,'--og','LineWidth', 3);
hold on ; grid on;
plot(t,Xs,'k','LineWidth', 3);
clear Xm Xs Xp Xsp stdXm stdXs t mR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);
Xp=reshape(LE(1:m*fr),fr,m);
Xsp=reshape(QE(1:m*fr),fr,m);
Xm=nanmean(Xp');
Xs=nanmean(Xsp');
t=[Datam(1,4):23, 0:Datam(1,4)-1];
figure(10)
subplot(1,1,1);
set(gca,'FontSize',10);
plot(t,Xm,'--or','LineWidth', 3);
hold on ; grid on;
plot(t,Xs,'y','LineWidth', 3);
legend('R_n Obs.','R_n Sim.','\lambdaE Obs.','\lambdaE Sim.')
xlabel('Hour'); ylabel('[W/m^2]')
clear Xm Xs Xp Xsp stdXm stdXs t mR


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nnd2=length(jDay);
T_H = sum(T_H,2); T_L=sum(T_L,2);
EIn_H = sum(EIn_H,2); EIn_L=sum(EIn_L,2);
EIn = EIn_H +EIn_L;
T_Hd= sum(reshape(T_H(1:nnd2*24),24,nnd2));
T_Ld= sum(reshape(T_L(1:nnd2*24),24,nnd2));
EGd= sum(reshape(EG(1:nnd2*24),24,nnd2));
EInd= sum(reshape(EIn(1:nnd2*24),24,nnd2));
QEd= mean(reshape(QE(1:nnd2*24),24,nnd2));


GPP_H =(NPP_H+RA_H)*Ccrown';
GPP_L =(NPP_L+RA_L)*Ccrown';
LAI_H2 = LAI_H*Ccrown';
LAI_L2 = LAI_L*Ccrown';
for i=1:365
    LE_doy(i)=mean(QEd(jDay==i));
    GPP_H_doy(i)=mean(GPP_H(jDay==i));
    GPP_L_doy(i)=mean(GPP_L(jDay==i));
    T_H_doy(i)=mean(T_Hd(jDay==i));
    T_L_doy(i)=mean(T_Ld(jDay==i));
    EG_doy(i)=mean(EGd(jDay==i));
    EIn_doy(i)=mean(EInd(jDay==i));
    LAI_H_doy(i)=mean(LAI_H2(jDay==i));
    LAI_L_doy(i)=mean(LAI_L2(jDay==i));
end
figure(11)
set(gcf,'color','w');
subplot(1,3,1)
plot(1:365,LAI_H_doy+LAI_L_doy,'-k','LineWidth',2)
hold on
plot(1:365,LAI_H_doy,'ob')
plot(1:365,LAI_L_doy,'xr')
xlabel('Doy')
ylabel('LAI [m^{2} m^{-2}]')
legend('Total','High','Low')
set(gca,'FontSize',12)
subplot(1,3,2)
plot(1:365,GPP_H_doy+GPP_L_doy,'-k','LineWidth',2)
hold on
plot(1:365,GPP_H_doy,'-b','LineWidth',1.5)
plot(1:365,GPP_L_doy,'-r','LineWidth',1.5)
xlabel('Doy')
ylabel('GPP [gC m^{-2} day{-1}]')
set(gca,'FontSize',12)
legend('Total','High','Low')
subplot(1,3,3)
plot(1:365,EG_doy+T_H_doy+T_L_doy+EIn_doy,'-k','LineWidth',2)
hold on
plot(1:365,T_H_doy,'-b','LineWidth',1.5)
plot(1:365,T_L_doy,'-r','LineWidth',1.5)
plot(1:365,EG_doy,'-m','LineWidth',1.5)
plot(1:365,EIn_doy,'-g','LineWidth',1.5)
xlabel('Doy')
ylabel('ET [mm day{-1}]')
legend('Total','High','Low','Ground','Interc.')
set(gca,'FontSize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Daily Cycles per month
anw_tit = 1;
didascalia={'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};
dt=1; mR =0;
for j=1:12
    %%%%%%%%%%%%%%%%%%%%%%%%%
    LEm=LE(Datam(:,2)==j);
    QEpm=QE(Datam(:,2)==j);
    Rnsm=Rns(Datam(:,2)==j);
    Rnpm=Rn(Datam(:,2)==j);
    Hsm=Hs(Datam(:,2)==j);
    Hpm=H(Datam(:,2)==j);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    n=length(Rnsm); fr=24/dt;
    %%%%%%%%%%%%%%%%%%%%
    m=floor(n/fr);
    Xp=reshape(Rnsm(1:m*fr),fr,m);
    Xsp=reshape(Rnpm(1:m*fr),fr,m);
    Xm=nanmean(Xp');
    Xs=nanmean(Xsp');
    stdXm=nanstd(Xp');
    stdXs=nanstd(Xsp');
    %t=1:fr;
    t=[Datam(1,4):23, 0:Datam(1,4)-1];
    figure(12)
    subplot(4,3,j);
    set(gca,'FontSize',10);
    plot(t,Xm,'-k','LineWidth', 1.5);
    hold on ; grid on;
    plot(t,Xs,'--k','LineWidth', 2);
    %plot(t,stdXm,'^r','LineWidth', 2);
    %plot(t,stdXs,'^b','LineWidth', 2);
    mR=max(mR,nanmax(nanmax(Xm,Xs))+50);
    clear Xm Xs Xp Xsp stdXm stdXs t
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m=floor(n/fr);
    Xp=reshape(LEm(1:m*fr),fr,m);
    Xsp=reshape(QEpm(1:m*fr),fr,m);
    Xm=nanmean(Xp');
    Xs=nanmean(Xsp');
    stdXm=nanstd(Xp');
    stdXs=nanstd(Xsp');
    %t=1:fr;
    t=[Datam(1,4):23, 0:Datam(1,4)-1];
    figure(12)
    subplot(4,3,j);
    set(gca,'FontSize',10);
    plot(t,Xm,'-r','LineWidth', 1.5);
    hold on ; grid on;
    plot(t,Xs,'--r','LineWidth', 2);
    %plot(t,stdXm,'^r','LineWidth', 2);
    %plot(t,stdXs,'^b','LineWidth', 2)
    clear Xm Xs Xp Xsp stdXm stdXs t
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m=floor(n/fr);
    Xp=reshape(Hsm(1:m*fr),fr,m);
    Xsp=reshape(Hpm(1:m*fr),fr,m);
    Xm=nanmean(Xp');
    Xs=nanmean(Xsp');
    stdXm=nanstd(Xp');
    stdXs=nanstd(Xsp');
    %t=1:fr;
    t=[Datam(1,4):23, 0:Datam(1,4)-1];
    figure(12)
    subplot(4,3,j);
    set(gca,'FontSize',10);
    plot(t,Xm,'-b','LineWidth', 1.5);
    hold on ; grid on;
    plot(t,Xs,'--b','LineWidth', 2);
    %plot(t,stdXm,'^r','LineWidth', 2);
    %plot(t,stdXs,'^b','LineWidth', 2);
    clear Xm Xs Xp Xsp stdXm stdXs t
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if anw_tit == 1
        if j <= 3
            title('Energy Budget')
        end
    else
        title(didascalia{j});
    end
    if j == 12
        legend('R_n OBS','R_n SIM','\lambda E OBS','\lambda E  SIM','H OBS','H SIM')
    end
    if j >=10
        xlabel('Hour');
    end
    if j == 1 || j== 4 || j == 7 || j == 10
        ylabel('[W/m^2]') ;
    end
    
end
for j=1:12
    figure(12)
    subplot(4,3,j);
    axis([0 24 -100 mR])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(111)
subplot(2,1,1)
plot(Date,GPP1,'--b')
hold on
plot(Date,GPP2,'--b')
plot(Date,GPP3,'--b')
plot(Date,GPP4,'--b')
plot(Date,(An_H+Rdark_H+An_L+Rdark_L)*Ccrown','k','LineWidth',1.5);
grid on
datetick('x',11)
title('GPP')
xlabel('Year'); ylabel('[umol m^{-2} s^{-1}]')
legend('OBS-1','OBS-2','OBS-3','OBS-4','SIM')



R=corrcoef(GPP(not(isnan(GPPs))),GPPs(not(isnan(GPPs))));

figure(112)
subplot(1,1,1)
set(gca,'Fontsize',11)
plot(GPP,GPPs,'sk','LineWidth',1,'Markersize',2)
hold on ; grid on;
title('GPP [umol m^{-2} s^{-1}]')
xlabel('SIM'); ylabel('OBS')
plot(0:60,0:60,'--y','LineWidth',1.5);
title(strcat('GPP R =', num2str(R(1,2))))



figure(113)
subplot(2,1,1)
plot(Dd,GPPd_m,'--r','LineWidth',1.5);
hold on ; grid on;
plot(Dd,GPPd_s,'b','LineWidth',1.5);
title('GPP')
datetick('x',11)
xlabel('Date'); ylabel('[gC m^{-2} day^{-1}]')
legend('OBS','SIM')
subplot(2,1,2)
plot(GPPd_s,GPPd_m,'sk','LineWidth',1.5);
hold on ; grid on;
title('GPP [gC m^{-2} day^{-1}]')
xlabel('SIM'); ylabel('OBS')
plot(0:15,0:15,'--y','LineWidth',1.5);
R=corrcoef(GPPd_s(not(isnan(GPPd_m))),GPPd_m(not(isnan(GPPd_m))));
title(strcat('GPP R =', num2str(R(1,2))))


yys=min(yyy)-1;

for yy=yyy
    GPPyr_m(yy-yys)=365*nanmean(GPPd_m(find(year(Dd)==yy)));
    GPPyr_s(yy-yys)=365*nanmean(GPPd_s(find(year(Dd)==yy)));
end

yy=yyy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(116)
subplot(1,2,1)
plot(yy,GPPyr_m,'--sr','LineWidth',1.5);
hold on ; grid on;
plot(yy,GPPyr_s,'--^b','LineWidth',1.5);
title('GPP')
xlabel('Year'); ylabel('[gC m^{-2} yr^{-1}]')
legend('OBS','SIM')
subplot(1,2,2)
plot(GPPyr_s,GPPyr_m,'sk','LineWidth',1.5);
hold on ; grid on;
title('GPP [gC m^{-2} yr^{-1}]')
xlabel('SIM'); ylabel('OBS')
plot(500:2500,500:2500,'--y','LineWidth',1.5);



if OPT_SoilBiogeochemistry == 1
    
    %%%%% ANNUAL GPP
    YY=year(Dd);
    for i=yyy
        I=year(Dd)==i;
        NEEym(i-yym1)=mean(NEEd(I));
        NEEys(i-yym1)=mean(NEEd_s(I));
    end
    
    figure(214)
    subplot(1,1,1)
    plot(yyy,NEEym(:),'-ob')
    a = polyfit(yyy,NEEym(:)',1);
    slm=a(1);
    hold on
    plot(yyy,NEEys(:),'-xr')
    a = polyfit(yyy,NEEys(:)',1);
    sls=a(1);
    xlim([min(yyy)-1 max(yyy)+1]);
    xlabel('Year'); ylabel('[gC m^{-2} day^{-1}]')
    title(strcat('OBS Trend=',num2str(slm),'; SIM Trend=',num2str(sls)))
    legend('OBS','SIM')
    
    
    %%% Trend for NEE - Seasonal
    for j=1:365
        for i=yyy
            I=year(Dd)==i &  jDay'==j ;
            NEE_doy(i-yym1)=mean(NEEd(I));
            NEE_doy2(i-yym1)=mean(NEEd_s(I));
        end
        M_NEE(j) = nanmean(NEE_doy);
        M_NEE_sim(j) = nanmean(NEE_doy2);
        clear  NEE_doy NEE_doy2
    end
    
    SLOPEp = M_NEE;
    window_S =30;
    SLO_F = 0*SLOPEp;
    V=[ SLOPEp , SLOPEp , SLOPEp];
    Vf= filter(ones(1,window_S )/window_S,1,V);
    Vf(1:window_S )=Vf(window_S +1);
    SLO_F=Vf([366:730]+window_S/2);
    M_NEE_smooth = SLO_F;
    %%%
    SLOPEp = M_NEE_sim;
    window_S =30;
    SLO_F = 0*SLOPEp;
    V=[ SLOPEp , SLOPEp , SLOPEp];
    Vf= filter(ones(1,window_S )/window_S,1,V);
    Vf(1:window_S )=Vf(window_S +1);
    SLO_F=Vf([366:730]+window_S/2);
    M_NEE_sim_smooth = SLO_F;
    
    
    figure(215)
    subplot(1,1,1)
    plot(1:365,M_NEE,'ob')
    hold on
    plot(1:365,M_NEE_sim,'xr')
    plot(1:365,M_NEE_smooth,'-b','Linewidth',2);
    plot(1:365,M_NEE_sim_smooth,'-r','Linewidth',2);
    xlabel('Doy'); ylabel('[gC m^{-2} day^{-1} year^{-1}]')
    legend('OBS','SIM')
    title('NEE')
    
    figure(211)
    subplot(2,1,1)
    plot(NEEd)
    hold on ; grid on ;
    plot(NEE,'r')
    title('NEE')
    legend('Mea','Sim')
    subplot(2,1,2)
    plot(REd)
    hold on ; grid on ;
    plot(RE,'r')
    title('RE')
    legend('Mea','Sim')
    
    figure(212)
    M=8;
    plot(NEE,NEEd,'x')
    axis([ -M M -M M])
    xlabel('Simulated')
    ylabel('Measured ')
    R=corrcoef(NEE(not(isnan(NEEd))),NEEd(not(isnan(NEEd))));
    hold on ;
    grid on ;
    x=-M:1:M;
    plot(x,x,'r')
    title(strcat('NEE R=', num2str(R(1,2))))
    
    r=0;
    for yy=yyy
        for k=1:12
            r=r+1;
            I= find((year(Dd)==yy)&(month(Dd)==k));
            NEE_m(r)=nanmean(NEEd(I));
            NEE_s(r)=nanmean(NEE(I));
            RE_m(r)=nanmean(REd(I));
            RE_s(r)=nanmean(RE(I));
        end
    end
    
    R=corrcoef(NEE_m,NEE_s); R2=R(1,2)^2;
    figure(213)
    subplot(2,1,1)
    plot(1:r,NEE_m,'--sr','LineWidth',1.5);
    hold on ; grid on;
    plot(1:r,NEE_s,'--^b','LineWidth',1.5);
    title('NEE')
    xlabel('Month'); ylabel('[gC m^{-2} day^{-1}]')
    legend('OBS','SIM')
    title(strcat('NEE R^2 =', num2str(R2)))
    R=corrcoef(RE_m,RE_s); R2=R(1,2)^2;
    subplot(2,1,2)
    plot(1:r,RE_m,'--sr','LineWidth',1.5);
    hold on ; grid on;
    plot(1:r,RE_s,'--^b','LineWidth',1.5);
    title('RE')
    xlabel('Month'); ylabel('[gC m^{-2} day^{-1}]')
    legend('OBS','SIM')
    title(strcat('RE R^2 =', num2str(R2)))
    
    clear NEE_m NEE_s  RE_m RE_s
    r=0;
    for yy=yyy
        r=r+1;
        I= find((year(Dd)==yy));
        NEE_m(r)=365*nanmean(NEEd(I));
        NEE_s(r)=365*nanmean(NEE(I));
        RE_m(r)=365*nanmean(REd(I));
        RE_s(r)=365*nanmean(RE(I));
    end
    R=corrcoef(NEE_m,NEE_s); R2=R(1,2)^2;
    
    yy=yyy;
    figure(216)
    subplot(2,1,1)
    plot(yy,NEE_m,'--sr','LineWidth',1.5);
    hold on ; grid on;
    plot(yy,NEE_s,'--^b','LineWidth',1.5);
    title('NEE')
    xlabel('Year'); ylabel('[gC m^{-2} day^{-1}]')
    legend('OBS','SIM')
    title(strcat('NEE R^2 =', num2str(R2)))
    R=corrcoef(RE_m,RE_s); R2=R(1,2)^2;
    subplot(2,1,2)
    plot(yy,RE_m,'--sr','LineWidth',1.5);
    hold on ; grid on;
    plot(yy,RE_s,'--^b','LineWidth',1.5);
    title('RE')
    xlabel('Year'); ylabel('[gC m^{-2} day^{-1}]')
    legend('OBS','SIM')
    title(strcat('RE R^2 =', num2str(R2)))
    
    
    
    
end



%%%%%%%%% SPECIFIC PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OPT_CASE == 3
    LAI_obs=[3.47455547298436,20100831
        3.44454656800764,20090831
        3.42459364277721,20080831
        3.28956981153964,20070831
        3.1217818069532,20060831
        3.00030703750454,20050831
        2.82268133773155,20040831
        2.61943344744443,20030831
        2.59373934594469,20020831
        3.29385279274002,20010831
        3.23747620551137,20000831
        3.17617142600597,19990831
        3.1348701438841,19980831
        3.10409988920112,19970831
        3.05070335087483,19960831
        3.01291225218796,19950831];
    Date_LAI=datenum(num2str(LAI_obs(:,2)),'yyyymmdd');
    LAI_obs=LAI_obs(:,1);
    
    Dd=Date(1):1:(Date(end)+1);
    
    figure(783)
    plot( Date_LAI,LAI_obs,'ko','markersize',2);
    hold on
    plot(Dd,(LAI_L+LAI_H)*Ccrown','--b','LineWidth',2.0);
    datetick('x','keepticks','keeplimits')
    title('LAI')
    ylabel('[m^2m^{-2}]')
    xlabel('Date')
    legend('OBS','SIM')
end
if OPT_CASE == 5
    load('E:\SIM_DISC\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Morgan_Monroe\Biometric_data_MMSF.mat')
    Dd=Date(1):1:Date(end)+1;
    figure(783)
    plot( Date_LAI_obs,LAI_obs,'ko','markersize',2);
    hold on
    plot(Dd,(LAI_L+LAI_H)*Ccrown','--m','LineWidth',2.0);
    datetick('x','keepticks','keeplimits')
    title('LAI')
    ylabel('[m^2m^{-2}]')
    xlabel('Date')
    legend('OBS','SIM')
end
if OPT_CASE == 9
    load('E:\SIM_DISC\ARCHIVE\Università deposito II\Fluxnet Project\Fluxnet - Data\Negrisia\Add_data_Negrisia.mat')
    Dd=Date(1):1:Date(end)+1;
    figure(784)
    plot( Date_LAI,LAI_obs,'ko','markersize',2);
    hold on
    plot(Dd,(LAI_L+LAI_H)*Ccrown','--m','LineWidth',2.0);
    datetick('x','keepticks','keeplimits')
    title('LAI')
    ylabel('[m^2m^{-2}]')
    xlabel('Date')
    legend('OBS','SIM')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for yy=yyy
        FP_yr(yy-yys)=365*nanmean(Sfr_H(find(year(Dd)==yy)));
        FP_oo(yy-yys)=sum(HR_crop(find(year(Date_HR_crop)==yy)));
    end
    yy=yyy;
    figure(786)
    subplot(1,1,1)
    plot(yy,FP_yr,'--sr','LineWidth',1.5);
    hold on
    plot(yy,FP_oo/2,'ko','LineWidth',1.5);
    hold on
    title('Harvest')
    xlabel('Year')
    ylabel('[gC m^{-2} yr^{-1}')
    legend('SIM','OBS')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %figure(787)
    %plot( Date_Ag_tot_prod,Ag_tot_prod,'ko','markersize',2);
    %hold on
    %plot(Date_HR_crop,HR_crop,'ko','markersize',2);
    %hold on
    %plot(Dd,(LAI_L+LAI_H)*Ccrown','--m','LineWidth',2.0);
    %datetick('x','keepticks','keeplimits')
    
end
%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



OSIM=zeros(length(Date),length(Idp));
for i=1:length(Date)
    OSIM(i,:)=interp1(Zs,O(i,:),Idp);
end
figure(301)
subplot(2,1,1)
plot(Date,SWC1/100,'xk','LineWidth',2.5);
hold on ; grid on ;
plot(Date,OSIM(:,1),'r','LineWidth',2.5);
datetick('x',11)
axis([ Date(1) Date(end) 0 0.8])
ylabel('\theta [-]')
title('a) Soil Moisture Depth # [cm]')
subplot(2,1,2)
plot(Date,SWC1/100,'xk','LineWidth',2.5);
hold on ; grid on ;
plot(Date,OSIM(:,2),'r','LineWidth',2.5);
datetick('x',11)
axis([ Date(1) Date(end) 0 0.8])
ylabel('\theta [-]')
title('b) Soil Moisture Depth # [cm]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(302)
%%%%%%%%%%%%
subplot(1,1,1)
plot(Date,Tdp,'r','LineWidth',1.0);
hold on ; grid on ;
plot(Date,Ts1,'--k','LineWidth',1.0);
hold on ; grid on ;
plot(Date,Ts2,'--m','LineWidth',1.0);
hold on ; grid on ;
title('Soil Temperature')
legend('OBS. # [cm]','OBS. # [cm]','SIM.')
xlabel('Date');
ylabel('[°C]')
axis([ Date(1) Date(end) min(Ts) max(Ts) ])
datetick('x',11)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(303)
TSIM=zeros(length(Date),length(Idp));
for i=1:length(Date)
    TSIM(i,:)=interp1(Zs,Tdp(i,:),Idp);
end
%%%%%%%%%%%%
subplot(2,1,1)
plot(Date,TSIM(:,1),'r','LineWidth',1.5);
hold on ; grid on ;
plot(Date,Ts1,'--k','LineWidth',1.0);
title('Soil Temperature # [cm]')
legend('SIM. ','OBS.')
xlabel('Date');
ylabel('[°C]')
axis([ Date(1) Date(end) min(Ts) max(Ts) ])
datetick('x',11)
subplot(2,1,2)
plot(Date,TSIM(:,2),'r','LineWidth',1.5);
hold on ; grid on ;
plot(Date,Ts2,'--k','LineWidth',1.0);
title('Soil Temperature  # [cm]')
legend('SIM. ','OBS.')
xlabel('Date');
ylabel('[°C]')
axis([ Date(1) Date(end) min(Ts) max(Ts) ])
datetick('x',11)


