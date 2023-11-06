% Plot RMSE statistics
% for forecast groups (i.e. initialized withe same hindcast)
%
% computed in calc_rmse_2Dssh.m
% or parallel version calc_rmse_2Dssh_paral.m
%
% RMSE between NEMO and 
% HYCOM analysis SSH fields
% Forecasts: use HYCOM hindcasts for initial fields
%  Forecasts: decision tree:
%
%  Hindcast group (which is used to initialize f/cst): 3, 7 or 8 currently
% 2 - not finished
% H/cast #2 - Full 2D SSH  -- not finished
%        #3 - AVISO SSH tracks only
%        #6 - GoM T/S profiles 30th points NEMO
%        #7 - AVISO + UGOS PIES (T/S profiles)
%        #8 - AVISO + extended PIES
%
%  Predictability experiments:
%        #10 - NEMO+GLORYS interpolated into HYCOM
%      |
%      V
%  TimePeriod during which the f/casts are initialized (May-June 2011, Jan-Feb 2012)
%     |
%     V
%  Forecasts (90 day f/csts shifted 7 days - 7 fcsts)
%
%


addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear



pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';

% Info about All hindcast experiments
load('hycom_tsis_expts.mat');

% Forecasts < 10 - OSSE f/casts
%                 each f/cast group: 7 f/casts with 1 week shift
%           10-16 - predictability with NEMO+GLORYS interpolated into HYCOM
%                   1 run initialized at day 0 then 1 week shifted
%                   1 f/cast = 1 run, each f/cast is 1 week shifted
%                   so 10-16 f/casts are similar to 1 f/cast group 
%                   This is 1a experiments
% 
%     Do not use 1b experiments - when 1a run is used +/- 1,2 days
%     to initialize f/casts with perturbed ICs around day 0 of 1a run
%     day 0 - starts May 8, then shifts 1 week for each 1a
%
%    So, here I use only 1a Predict runs 10-16 as 1 f/cast group
%

%Z0 = -10; % discard close to coastline points
Z0 = -200;  % discard shelf to reduce near-coastal noise

Ntime = 2;
IFCST=[2, 3, 6, 7, 8, 10]; % forecast groups
Nfgr=length(IFCST); % # of forecasts groups 

%
% Combine RMSE fields
% for OSSE f/casts
RMSE = sub_combine_RMSE_OSSEfcst(IFCST,pthout,EXPT,Z0);
%
% Combine Predict1a experiments - f/cast #10 only
RMSE = sub_combine_RMSE_prdct1a(RMSE,pthout,Z0);

% Weekly pools
WK=[1:7:92]; % weekly pools
nwk=length(WK)-1;

% Pool all runs for the same f/cast group
% Average by months:

for itime=1:2
  for ifc=1:Nfgr
    POOL(ifc).Time(itime).pm1=[];
    POOL(ifc).Time(itime).pm2=[];
    POOL(ifc).Time(itime).pm3=[];
    PRST(ifc).Time(itime).pm1=[];
    PRST(ifc).Time(itime).pm2=[];
    PRST(ifc).Time(itime).pm3=[];



    for iwk=1:nwk
      POOL(ifc).Time(itime).WK(iwk).pw=[];
      PRST(ifc).Time(itime).WK(iwk).pw=[];
    end

  end
end

for itime=1:Ntime
  for ifc=1:Nfgr
%
% Pool all runs - spatial mean RMSE
% groups by months
    dmm = RMSE(ifc).Time(itime).RMSE_mean;
    POOL = sub_pool_month_week(POOL,dmm,ifc,itime,WK);

% Persistence:
    dmm = RMSE(ifc).Time(itime).RMSEprst_mean;
    PRST = sub_pool_month_week(PRST,dmm,ifc,itime,WK);
  end
end

% Hindcast colors 
% match with hindcast for
% comparison
CLR = [0.4 0.8 0.2;...
       0   0.6 0.9; ...
       0.5   1   0.7; ...
       1   0.4 0.5; ... 
       0.  1  0.3; ...
       1   0   0.9; ...
       0.8 0   0.4; ...
       1   0.8 0; ...
       0.8 0.5 0; ...
       0.7 0.6 0.4; ...
       0.5 0.2 1];


btx='plot_RMSEbar_OSSEfcst_Prdct.m';


% Plot: all F/csts in 1 group for 3 30-day time periods
% coordinates of the Bars relative to X=# of the 30-day period, 1,2,3
nfg=40;
anls_nm=sprintf('RMSE(m) SSH months');
%sub_plotMHD_barsrt_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
%sub_plotMHDprst_barsrt_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,PRST,anls_nm);
% Plot persistence as separate bar next to expt:
%sub_plotMHD_prstbar_barsrt_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,PRST,anls_nm);
% Plot perst as a line IQR with median dot
sub_plotMHD_prstline_barsrt_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,PRST,anls_nm);

bottom_text(btx,'pwd',1);

% Plot weekly bars
f_week=0;
if f_week==1
  nfg=41;
  anls_nm=sprintf('RMSE(m) SSH weeks');
%sub_plotMHD_barsrt_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,WK,anls_nm);
   sub_plotMHDprst_barsrt_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,PRST,WK,anls_nm);
end

bottom_text(btx,'pwd',1);








