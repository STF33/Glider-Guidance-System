% Example how to plot RMSE 
% of T anomalies at 200 m 
% from the HYCOM predicatbility experiments 1b
% 
%  Predictability experiments:
%  NEMO+GLORYS interpolated into HYCOM
%      #10 - 01/05/2011 (Time 1) and 01/01/2012 (Time 2)
%      #11 - 08/05/2011 (Time 1) and 08/01/2012 (Time 2)
% etc
% Predictability experiments: iFcst 10, 11, ..., 16:
% Forecasts: use interpolated fields as initial fields
%            within each f/cast - perorm additional f/cast runs
%            initialized from day t0=7 days of HYCOM main f/cast
%            run1 = main f/cast, 
%            run2 = -2day shift used as IC at day t0
%            run3 = -1 day shift, etc.
%            run4 = +1 days
%            run5 = +2 days

% Select experiment iFcst = 10, ..., 16
% perturbation run: irun = 2, ..., 5
% Time period: itime = 1, 2 - 2011/2012

Z0    = -200;
itime = 1; % Runs in 2011
iFcst = 10;  % The run that starts t0= May 08 2011
irun  = 3;    % IC is -1 day shift wrt to t0

pthout = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
frmseout = sprintf('%sRMSE_Prdct1b_dT%4.4i_%s.mat',pthout,abs(Z0),fcstname);
load(frmseout);

% Pick day to plot RMSE
itime=30; % day 30 of the f/cast

% Note that saved domain is subdomain of the whole HYCOM computational domain
HH   = RMSERR.HH_gom; % topography
Iocn = RMSERR.Indx_gom; % linear indices of T points insed GoM
RMSE = sqrt(RMSERR.ERR_squared);  % RMSE
rmse = RMSE(:,itime);

A = HH*nan;
A(Iocn) = rmse;
figure(1); clf;
pcolor(A); shading flat;
colorbar
caxis([0 2.5]);



