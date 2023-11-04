% Plot RMSE statistics
% for OSE 3mo forecasts
% RMSE
% RMSE computed for OSE vs 50.1 GOMu reanalysis
% RMSE computed at calc_rmse_OSEhindcast_501GOMu.m
% Annual mean
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

btx = 'plot_RMSEbar_OSEhindcast.m';

% for 2009 and 2011 - not full year - only covered by PIES obs
clear RMSE
RMSE = sub_combine_RMSE_OSEhndcst(pthout,Z0);


% Hindcast colors 
% match with hindcast for
% comparison
CLR = [0,0.3,0.7; ...
       0.7,0.4,0];

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

  for iyr=1:3
    mdn = RMSE(ifc).median(iyr);
    p25 = RMSE(ifc).p25(iyr);
    p75 = RMSE(ifc).p75(iyr);

    ixx = iyr;
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

xlabel('Years');

set(gca,'tickdir','out',...
        'xlim',[0.5 3.5],...
        'xtick',[1:3],...
        'xticklabel',{'2009','2010','2011'},...
        'ylim',[0 0.25],...
        'ytick',[0:0.05:0.3],...
        'ygrid','on',...
        'Fontsize',12);

title('Median & IQR SSH RMSE (m) OSE analysis');



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





