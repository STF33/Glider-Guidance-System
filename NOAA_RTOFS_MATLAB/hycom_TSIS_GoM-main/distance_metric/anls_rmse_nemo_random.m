% Analyze RMSE statistics for random RMSE
% see calc_rmse_ssh_OSSErandom.m

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

rmse = [];
for YR=2011:2012
  frmseout = sprintf('%sRMSE_NEMO_random_%i.mat',pthout,YR);
  fprintf('Loading %s\n',frmseout);
  load(frmseout);

  dmm  = RMSERR.ERR_squared;
  amm  = sqrt(mean(dmm,1));
  amm  = amm(:);
  rmse = [rmse;amm];
end


dx = 0.005;
Xf = [0.14:dx:0.24];
Nf = hist(rmse,Xf);
Nf = Nf/sum(Nf);

figure(1); clf;
axes('Position',[0.1,0.5,0.8,0.4]);
hb = bar(Xf,Nf,0.98);
set(hb,'Facecolor',[0.6, 0.6, 0.6]);
hold on;

% plot IQR
xmd = median(rmse);
xlw = prctile(rmse,25);
xup = prctile(rmse,75);
dbr = 0.01;
ybr = 0.14;

plot(xmd,ybr,'.','Color',[0 0 0],'Markersize',25);
plot([xlw, xup],[ybr,ybr],'k-','linewidth',2);
plot([xlw,xlw],[ybr-dbr, ybr+dbr],'k-','linewidth',2);
plot([xup,xup],[ybr-dbr, ybr+dbr],'k-','linewidth',2);

set(gca,'tickdir','out',...
        'ytick',[0:0.05:0.5],...
        'xtick',Xf,...
        'xlim',[0.15 0.23],...
        'ylim',[0 0.17],...
        'ygrid','on');

title('RMSE from random NEMO NR SSH 2011-2012');

axes('Position',[0.1,0.1,0.8,0.35]);
stxt = sprintf('IQR = %5.3f/%5.3f, median = %5.3f, mean =%5.3f',...
                xlw,xup,xmd,mean(rmse));

text(0.1,0.1,stxt);
set(gca,'visible','off');

btx = 'anls_rmse_nemo_random.m';
bottom_text(btx,'pwd',1);


