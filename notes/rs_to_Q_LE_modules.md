It seems that in T&C, Q_LE is called QE.

# Path
MAIN_FRAME.m -> HYDROLOGIC_UNIT.m -> SVAT_UNIT -> Heat_fluxes

Heat_fluxes computes Q_LE (as QE) similarly to what I did in my python pipeline, but with scaled terms (by L, Ls?)
It uses terms T_H, T_L that depend on Ccrown, and Tpot_{H, L} (and are then passed through a min).
Tpot_{L,H} depend on Tpot_{L,H}{sun,shd} which is expressed as a complicated function of Cice, Csno, ro, qTvLsun_sat, qTa, ra, rb_L, LAI_L, FsunL, dw_L, rap_H, and finally rs_sunL.
The question now, is whether any of those predictors depend on rs. If that's not the case, I could easily implement this part and fetch the predictors from the T&C output. 