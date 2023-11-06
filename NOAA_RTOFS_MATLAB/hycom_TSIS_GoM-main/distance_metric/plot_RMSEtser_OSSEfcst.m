% Plot RMSE time series for 3 months OSE/OSSE f/casts with persistence
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
RMSE(end).Hnd_name = 'NEMO_Intrp';

% 
% Combine RMSE from OSE PIES/noPIES f/casts:
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
Z0 = -200; % only deep GoM
RMSE_ose = sub_combine_RMSE_OSEfcst(pthout,Z0);

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

%clrF = [0,0,0];  % OSSE f/cast color
% PIES/noPIES OSEs:
CLRPI = [0,0.4,1; ...
       1,0.6,0];
clrpi   = CLRPI(1,:);
clrnopi = CLRPI(2,:); 

btx='plot_RMSEtser_OSSEfcst.m';

POS = [0.1,0.6,0.7,0.3; ...
       0.1,0.15, 0.7,0.3];
  
Rsat = 0.178;
%keyboard
for ik = 1:Nfgr+1
  figure(ik); clf;
  set(gcf,'Position',[1647, 720, 800, 600]);
  
  for itime=1:2
    pos = POS(itime,:);
    axes('Position',pos);
    hold on;
    if ik<=Nfgr
      hnm = RMSE(ik).Hnd_name;
      drr = RMSE(ik).Time(itime).RMSE_mean;
      dmm = RMSE(ik).Time(itime).RMSEprst_mean;  % persist
    else
% OSE
      hnm = RMSE_ose(itime).Hnd_name;
      drr = RMSE_ose(itime).RMSE_mean;
      dmm = RMSE_ose(itime).RMSEprst_mean;
    end
    rr  = mean(drr,2);
    [J,I] = find(drr >= Rsat); 
    prdct_rr = min(J);         % day predictability
    if isempty(J), prdct_rr = 99999; end; 
    rrp = mean(dmm,2); 
    [J,I] = find(dmm >= Rsat); 
    prdct_prst = min(J);         % day predictability
    if isempty(J), prdct_prst = 99999; end;

%    clrF = CLRPI(itime,:);
    clrF = [0 0.2 1];
    clrP = [0 0 0];
    plot(drr,'Linewidth',1,'Color',[0.5 0.6 1]);
    plot(dmm,'Linewidth',1,'Color',[0.7 0.7 0.7]);
    plot(rr,'Linewidth',2,'Color',clrF);        % f/cast
    plot(rrp,'-','Linewidth',2,'Color',clrP);  % persist
%    plot(rrpi,'Linewidth',2,'Color',clrpi);
%    plot(rrpi_ps,'--','Linewidth',2,'Color',clrpi);
%    plot(rrnopi,'Linewidth',2,'Color',clrnopi);
%    plot(rrnopi_ps,'--','Linewidth',2,'Color',clrnopi);

% Saturation value:
    plot([0 91],[Rsat Rsat],'Linewidth',1.2,'Color',[0.9 0.2 0]);

    set(gca,'tickdir','out',...
            'ylim',[0 0.25],...
            'xlim',[1 90],...
            'ytick',[0:0.05:0.3],...
            'xtick',[0:10:100],...
            'xgrid','on',...
            'ygrid','on');

    xlabel('F/cast days');
    iyr = 2011;
    if itime==2, iyr=2012; end;
    stl = sprintf('RMSE f/cast %s %i, Dpred=%i, Dpred_prst=%i',hnm,iyr,prdct_rr,prdct_prst);
    title(stl,'Interpreter','none');
  end

  bottom_text(btx,'pwd',1);

  % Legend:
  axes('Position',[0.85 ,0.55,0.12,0.3]);
  hold on;
  x1=0;
  x2=0.2;
  y1=1;
  plot([x1,x2],[y1,y1],'-','Color',clrF,'linewidth',2);
  text(x2+0.1,y1,'E[F/cast]','Fontsize',12);

  y1=2;
  plot([x1,x2],[y1,y1],'-','Color',clrP,'linewidth',2);
  text(x2+0.1,y1,'E[Persist]','Fontsize',12);

  set(gca,'xlim',[0 0.8],...
          'ylim',[0. 3.],...
          'visible','off');
end









