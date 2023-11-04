% Plot MHD time series for OSE/OSSE forecasts
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


pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';

%IFCST=[2, 3, 6, 7, 8, 10]; % forecast groups
IFCST=[2];
Nfgr=length(IFCST); % # of forecasts groups 

pthmat = pthmat1;

% All hindcast experiments
load('hycom_tsis_expts.mat');

%
% Combine MHD by forecast groups
% for OSSE f/casts
ifc=0;
MHD = sub_combine_MHD_OSSEfcst(IFCST,pthmat);

% Combine all Predictability 1a experiments into 1 f/cast group #10
MHD = sub_combine_MHD_prdct1a(MHD,pthmat1);
MHD(end).Time(1).Fcst_name = 'NEMO_Intrp';
MHD(end).Time(2).Fcst_name = 'NEMO_Intrp';

% Combine MHD from OSE PIES/noPIES f/casts:
AA = sub_combine_MHD_OSEfcst; 

% Reshape struct similar to OSSE, where time 1 = PIES, time 2 = noPIES
a1 = [];
a2 = [];
MHD_ose = struct;
MHD_ose.Time(1).Fcst_name = 'PIES';
MHD_ose.Time(2).Fcst_name = 'noPIES';
l1 = 91;
c1 = 0;
c2 = 0;
c3 = 0; 
c4 = 0;
nll = length(AA);
for ill=1:nll
  esim = AA(ill).esim;
  dmm  = AA(ill).mhd;
  if length(dmm) < l1
    dmm(end+1:l1) = nan;
  end
  switch(esim)
   case('PIES');
    c1 = c1+1;
    MHD_ose.Time(1).MHD(1:l1,c1) = dmm(1:l1);
   case('Persist PIES');
    c2 = c2+1;
    MHD_ose.Time(1).MHDprst(1:l1,c2) = dmm(1:l1);
   case('noPIES');
    c3 = c3+1;
    MHD_ose.Time(2).MHD(1:l1,c3) = dmm(1:l1);
   case('Persist noPIES');
    c4 = c4+1;
    MHD_ose.Time(2).MHDprst(1:l1,c4) = dmm(1:l1);
  end
end
  

% PIES/noPIES OSEs:
CLRPI = [0,0.4,1; ...
       1,0.6,0];
clrpi   = CLRPI(1,:);
clrnopi = CLRPI(2,:); 

btx='plot_MHDtser_OSSEfcst.m';

POS = [0.1,0.6,0.7,0.3; ...
       0.1,0.15, 0.7,0.3];

  
Rsat = 66.3;
nl = 0;
for ik = 1:Nfgr+1
  figure(ik); clf;
  set(gcf,'Position',[1647, 720, 800, 600]);
  
  for itime=2:2
    pos = POS(itime,:);
    axes('Position',pos);
    hold on;
    if ik<=Nfgr
      hnm = MHD(ik).Time(itime).Fcst_name;
      drr = MHD(ik).Time(itime).MHD;
      dmm = MHD(ik).Time(itime).MHD_prst;
    else
% OSE
      hnm = MHD_ose.Time(itime).Fcst_name;
      drr = MHD_ose.Time(itime).MHD;
      dmm = MHD_ose.Time(itime).MHDprst;
    end
% Smooth MHD time series
% Butterworth filter
    if length(drr) ~= nl
      nl  = length(drr);
      iid = 20;    % days for filtering
      Wn  = iid/nl;
      [Bf,Af] = butter(4,Wn,'low');
    end

    [a1,a2] = size(drr);
    d0 = nanmean(nanmean(drr));
    drr(isnan(drr)) = d0;
    d0 = nanmean(nanmean(dmm));
    dmm(isnan(dmm)) = d0;
    for kk=1:a2
      aa = drr(:,kk);
      af = filtfilt(Bf,Af,aa);
      drr(:,kk) = af;

      bb = dmm(:,kk);
      bf = filtfilt(Bf,Af,bb);
      dmm(:,kk) = bf;
    end

    rr  = nanmean(drr,2);
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

% Saturation value:
    plot([0 91],[Rsat Rsat],'Linewidth',1.2,'Color',[0.9 0.2 0]);

    set(gca,'tickdir','out',...
            'ylim',[0 150],...
            'xlim',[1 90],...
            'ytick',[0:25:200],...
            'xtick',[0:10:100],...
            'xgrid','on',...
            'ygrid','on');

    xlabel('F/cast days');
    iyr = 2011;
    if itime==2, iyr=2012; end;
    stl = sprintf('MHD f/cast %s %i, Dpred=%i, Dpred_prst=%i',hnm,iyr,prdct_rr,prdct_prst);
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









