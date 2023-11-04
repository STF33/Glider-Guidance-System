!/bin/csh
#
## Edit: input archive setenv D
# dates of the created restart
#
set echo
# --- Form a HYCOM restart file from a HYCOM archive file.
# # --- output is HYCOM restart files
# #
# # --- R is region name.
# # --- E is expt number.
# # --- X is decimal expt number.
# # --- T is topography number.
# # --- Y is year
# #
#
setenv OS Linux
setenv R IASx0.03
setenv DS  /gpfs/research/coaps/home/ddmitry/HYCOM_TSIS/IASx0.03
setenv S  ${DS}/data/restart
setenv SRC /nexsan/people/ddmitry/Net_ocean/HYCOM/hycom/ALL4
cd ${S}/.


setenv Y 001
setenv IDM    1401
setenv JDM    891
setenv KDM    41
setenv YRFLAG 3
setenv THBASE 34
setenv BACLIN 180




