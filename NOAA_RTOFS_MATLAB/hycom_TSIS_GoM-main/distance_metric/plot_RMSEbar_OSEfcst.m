% Plot RMSE statistics
% for OSE 3mo forecasts
% RMSE
% old: computed in calc_rmse_ssh_OSEfcst.m vs OSE hindcast SSH
% new: computed in calc_rmse_ssh_OSEfcst_501GOMu.m vs 50.1 GOMl reanalysis
%
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

esim1 = 'PIES';
esim2 = 'noPIES';

%
% Specify depth limit for RMSE calculation
Z0 = -200; % only deep GoM

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';

btx = 'plot_RMSEbar_OSEfcst.m';

%
% Combine RMSE fields
% for OSSE f/casts
RMSE = sub_combine_RMSE_OSEfcst(pthout,Z0);
%

% Weekly pools
WK=[1:7:92]; % weekly pools
nwk=length(WK)-1;

% Pool all runs for the same f/cast group
% Average by months:
POOL = struct;
for ifc=1:2
  for imo=1:3
    POOL(ifc).mo(imo).pm=[];
    PRST(ifc).mo(imo).pm=[];
  end

  for iwk=1:nwk
    POOL(ifc).WK(iwk).pw=[];
    PRST(ifc).WK(iwk).pw=[];
  end

end

for ifc=1:2
%
% Pool all runs - spatial mean RMSE
% groups by months
  dmm = RMSE(ifc).RMSE_mean;
  POOL = sub_pool_month_week_OSE(POOL,dmm,ifc,WK);

% Persistence:
  dmm = RMSE(ifc).RMSEprst_mean;
  PRST = sub_pool_month_week_OSE(PRST,dmm,ifc,WK);
end

% Hindcast colors 
% match with hindcast for
% comparison
CLR = [0,0.3,0.7; ...
       0.7,0.4,0];
clrp = [0.5,0.5,0.5];  % persistence
btx='plot_RMSEbar_OSEfcst.m';


% -------------------------
% Only monthly bars
% skip weekly
% ---------------------
for ifc=1:2
  for imo=1:3
    dmm = POOL(ifc).mo(imo).pm;
    mdn = nanmedian(dmm);
    p25 = prctile(dmm,25);
    p75 = prctile(dmm,75);

    POOL(ifc).mo(imo).Mdn = mdn;
    POOL(ifc).mo(imo).p25 = p25;
    POOL(ifc).mo(imo).p75 = p75;
  end

  for imo=1:3
    dmm = PRST(ifc).mo(imo).pm;
    mdn = nanmedian(dmm);
    p25 = prctile(dmm,25);
    p75 = prctile(dmm,75);

    PRST(ifc).mo(imo).Mdn = mdn;
    PRST(ifc).mo(imo).p25 = p25;
    PRST(ifc).mo(imo).p75 = p75;
  end
end


figure(1); clf;
set(gcf,'Position',[1573, 560, 949, 751]);
axes('Position',[0.1 0.5 0.6 0.3]);
hold on;

dx0 = 0.45;
for ifc=1:2
  clr=CLR(ifc,:);
  if ifc==1
    dx = -dx0;
  else
    dx = dx0;
  end

  for imo=1:3
    mdn = POOL(ifc).mo(imo).Mdn;
    p25 = POOL(ifc).mo(imo).p25;
    p75 = POOL(ifc).mo(imo).p75;

    ixx = imo;
    patch([ixx+dx ixx+dx ixx ixx],[0 mdn mdn 0],clr,'Edgecolor','none');

% Range:
    llw = p25;
    lup = p75;
    ix0=ixx+dx/2;
    plot([ix0-0.05 ix0+0.05],[llw llw],'k-','linewidth',2);
    plot([ix0-0.05 ix0+0.05],[lup lup],'k-','linewidth',2);
    plot([ix0 ix0],[llw lup],'k-','linewidth',1.6);
  end
end

% Plot persistence:
%CLRP = [0.,0.6,1; ...
%       1,0.6, 0];
clrP = [0.7 0.7 0.7];
for ifc=1:2
%  clr=CLRP(ifc,:);
  if ifc==1
    dx = -dx0/2;
  else
    dx = dx0*3/2;
  end

  for imo=1:3
    mdn = PRST(ifc).mo(imo).Mdn;
    p25 = PRST(ifc).mo(imo).p25;
    p75 = PRST(ifc).mo(imo).p75;
    dbr = 0.02;

    ixx=imo;
%    ix0=ixx+dx/2;
    ixM=ixx+dx/2;
% Persistence:
    plot(ixM,mdn,'.','Color',clrP,'Markersize',20);
    plot([ixM-dbr ixM+dbr],[p25 p25],'-','Color',clrP,'linewidth',1.5);
    plot([ixM-dbr ixM+dbr],[p75 p75],'-','Color',clrP,'linewidth',1.5);
    plot([ixM ixM],[p25 p75],'-','Color',clrP,'linewidth',1.5);
%    plot([ixx+dx ixx],[mdn mdn],'-','Color',clr,'Linewidth',3);
%    plot([ixx+dx ixx],[p25 p25],'-','Color',clr,'Linewidth',1.6);
%    plot([ixx+dx ixx],[p75 p75],'-','Color',clr,'Linewidth',1.6);
%    plot([ix0 ix0],[p25 p75],':','Color',clr);
  end
end
xlabel('Forecast Months');

set(gca,'tickdir','out',...
        'xlim',[0.5 3.5],...
        'xtick',[1:3],...
        'ygrid','on',...
        'Fontsize',12);

title('Median & IQR SSH RMSE (m) OSE 3-mo forecasts and Persistence');



bottom_text(btx,'pwd',1);

% Legend:
axes('Position',[0.1,0.1,0.4,0.25]);
hold on;
x1=0;
x2=1;
y1=1;
clr=CLR(1,:);
plot([x1,x2],[y1,y1],'-','Color',clr,'linewidth',20);
text(x2+0.1,y1,'PIES','Fontsize',14);

y1=2;
clr=CLR(2,:);
plot([x1,x2],[y1,y1],'-','Color',clr,'linewidth',20);
text(x2+0.1,y1,'noPIES','Fontsize',14);

set(gca,'xlim',[0 1.5],...
        'ylim',[0.8 2.3],...
        'visible','off');





