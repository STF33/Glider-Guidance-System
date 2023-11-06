% Example how to plot RMSE 
% of T anomalies at 200 m 
% from the HYCOM-NEMO predicatbility experiments 1a
% 
%  Predictability experiments:
%  NEMO+GLORYS interpolated into HYCOM
%      #10 - 01/05/2011 (Time 1) and 01/01/2012 (Time 2)
%      #11 - 08/05/2011 (Time 1) and 08/01/2012 (Time 2)
% etc
% Predictability experiments: iFcst 10, 11, ..., 16:
% THere is only 1 free run in these predictability experiments, irun=1

% Select experiment iFcst = 10, ..., 16
% irun = 1
% Time period: itime = 1, 2 - 2011/2012

Z0    = -200;
itime = 1; % Runs in 2011
iFcst = 10;  % The run that starts t0= May 08 2011
irun  = 1;   

pthout = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
frmseout = sprintf('%sRMSE_Prdct_dT%4.4i_%s.mat',pthout,abs(Z0),fcstname);
load(frmseout);

% Pick day to plot RMSE
itime=15; % day 30 of the f/cast

% Load HYCOM grid and topo:
pthdata = '/Net/kronos/ddmitry/hycom/TSIS/data2ncsu/';
fgrd = sprintf('%sHYCOM_grid_Iocn0200m.mat',pthdata);
load(fgrd);
HH   = GRID.topo;
Iocn = GRID.Iocn;
LON  = GRID.LON;
LAT  = GRID.LAT;

% Note that saved domain is subdomain of the whole HYCOM computational domain
RMSE = sqrt(RMSERR.ERR_squared);  % RMSE
rmse = RMSE(:,itime);

A = HH*nan;
A(Iocn) = rmse;
figure(1); clf;
pcolor(LON,LAT,A); shading flat;
hold on;
contour(LON,LAT,HH,[0 0],'k');
caxis([0 2.5]);
axis('equal');
set(gca,'xlim',[-98 -82],...
        'ylim',[18 30]);
colorbar


