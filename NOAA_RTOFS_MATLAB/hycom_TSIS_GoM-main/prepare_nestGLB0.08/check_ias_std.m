% Check tops_ias_std.nc file
% created for tsis assimilation
% see: /home/ddmitry/codes/PYTHON/hycom_TSIS/write_ias_std.py
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

pthin = '/nexsan/people/ddmitry/hycom/TSIS/IASx0.03/tsis_data/';
fnc1 = sprintf('%stops_ias_std_41l.nc',pthin);

pin = nc_varget(fnc1,'pin');
kk = 20;
AA = squeeze(pin(kk,:,:));

figure(1); clf;
pcolor(AA); shading flat;




