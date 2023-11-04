% Plot RMSE maps for OSE hindcasts
% with PIES vs no PIES (real observations)
% RMSE computed for OSE vs 50.1 GOMu reanalysis
% RMSE computed at calc_rmse_OSEhindcast_501GOMu.m
% Annual mean

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

yr = 2011;
%esim1 = 'PIES';
esim1 = 'noPIES';
esim2 = '501GOMu';
esimH = esim1;


%
% Specify depth limit for RMSE calculation
Z0 = -200; % only deep GoM

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';



% HYCOM:
rg=9806;  % convert pressure to depth, m
huge=1e20;

%Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
HH(isnan(HH))=100;


% GoM region HYCOM:
GOM=[366   489
   476   531
   583   560
   576   646
   508   827
   336   848
   204   829
    64   798
    19   746
    16   662
    12   578
    25   455
    71   382
   165   356
   281   400];

[XM,YM]=meshgrid([1:n],[1:m]);
INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
Iocn = find(INH==1 & HH<Z0);
Nocn=length(Iocn);
clear XM YM


% 
% Subsample to a smaller domain:
xnas1 = min(LON(Iocn));
xnas2 = max(LON(Iocn));
ynas1 = min(LAT(Iocn));
ynas2 = max(LAT(Iocn));
[J,I]=find(LON<xnas1);
its1=max(I);
[J,I]=find(LON>xnas2);
its2=min(I);
[J,I]=find(LAT<ynas1);
jts1=max(J);
[J,I]=find(LAT>ynas2);
jts2=min(J);

xt1=LON(jts1,its1);
yt1=LAT(jts1,its1);
xt2=LON(jts2,its2);
yt2=LAT(jts2,its2);


%frmseout = sprintf('%sRMSE%4.4im_OSEhindcastGOMu_%i.mat',pthout,abs(Z0),yr);
frmseout = sprintf('%sRMSE%4.4im_OSEhindcast%s_GOMu_%i.mat',pthout,abs(Z0),esimH,yr);

clear RMSERR
fprintf('Loading %s\n',frmseout);
load(frmseout);

dmm=sqrt(RMSERR.ERR_squared);

% Mean over PIES time period:
if yr == 2009
  RMSE = mean(dmm(:,91:end),2);  % Apr - Dec
elseif yr == 2011
  RMSE = mean(dmm(:,1:334),2);  % Jan - Nov.
else
  RMSE = mean(dmm,2);
end
btx = 'plot_RMSE_OSEhindcast_501GOMu_map.m';

stl = sprintf('RMSE SSH, OSEs %s vs 50.1GOMu rnls, Z0=%im,  %i',esim1,Z0,yr);

nfg=1;
sub_plot_rmse(nfg,HH,LON,LAT,Iocn,RMSE,stl);
bottom_text(btx,'pwd',1);


