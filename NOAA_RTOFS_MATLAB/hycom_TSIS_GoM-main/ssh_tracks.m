% SSH tracks
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

pth='/Net/kronos/ddmitry/hycom/TSIS/qcobs_WITHOUT_PITES/';
fnc = sprintf('%stsis_obs_ias_2009122700.nc',pth);

ssh = nc_varget(fnc,'ssh');
ssha = nc_varget(fnc,'av_ssh');



