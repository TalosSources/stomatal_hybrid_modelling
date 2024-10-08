##############################################################
### CLOSURE CO2 Concentration inside the stomatal           #
##############################################################
##############################
#function[DCi]=CO2_Concentration(Ci,IPAR,Ca,ra,rb,Ts,Pre,Ds,...
#    O,Owp,Oss,CT,VRmax,DS,Ha,FI,Oa,Do,a1,go,gmes,rjv)
###
#[CiF]= PHOTOSYNTESIS(Ci,IPAR,Ca,ra,rb,Ts,Ts,Pre,Ds,...
#    O,Owp,Oss,...
#    CT,VRmax,NaN,NaN,DS,Ha,FI,Oa,Do,a1,go);
#DCi = Ci - CiF;
#end

# Actually, it's a simple wrapper of photosynthesis_biochemical
###################################

from photosynthesis_biochemical import photosynthesis_biochemical

def CO2_Concentration(Cc,IPAR,Csl,ra,rb,Ts,Pre,Ds,Psi_L,Psi_sto_50,Psi_sto_99,CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,gmes,rjv):
     CcF,An,rs,Rdark,F755nm,GAM,gsCO2 = photosynthesis_biochemical(Cc,IPAR,Csl,ra,rb,Ts,Pre,Ds,Psi_L,Psi_sto_50,Psi_sto_99,CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,gmes,rjv)
     DCi = Cc - CcF
     return DCi