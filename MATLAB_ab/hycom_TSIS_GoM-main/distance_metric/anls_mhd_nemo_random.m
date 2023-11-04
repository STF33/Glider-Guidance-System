% Analyze RMSE statistics for random RMSE
% see calc_rmse_ssh_OSSErandom.m

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';

MHD = [];
for YR = 2011:2012
  fmat1 = sprintf('%sMHD_LCLCE_nemo_random_%i.mat',pthmat,YR);
  fprintf('Loading %s\n',fmat1);
  A = load(fmat1);
  MHD = [MHD;A.MHD];
end

dx = 5;
Xf = [45:dx:105];
Nf = hist(MHD,Xf);
Nf = Nf/sum(Nf);

figure(1); clf;
axes('Position',[0.1,0.5,0.8,0.4]);
hb = bar(Xf,Nf,0.98);
set(hb,'Facecolor',[0.6, 0.6, 0.6]);
hold on;

% plot IQR
xmd = median(MHD);
xlw = prctile(MHD,25);
xup = prctile(MHD,75);
dbr = 0.015;
ybr = 0.27;

plot(xmd,ybr,'.','Color',[0 0 0],'Markersize',25);
plot([xlw, xup],[ybr,ybr],'k-','linewidth',2);
plot([xlw,xlw],[ybr-dbr, ybr+dbr],'k-','linewidth',2);
plot([xup,xup],[ybr-dbr, ybr+dbr],'k-','linewidth',2);

set(gca,'tickdir','out',...
        'ytick',[0:0.05:0.5],...
        'xtick',Xf,...
        'xlim',[42 107],...
        'ylim',[0 0.30],...
        'ygrid','on');

title('MHD from random NEMO NR LC/LCE contours');

axes('Position',[0.1,0.1,0.8,0.35]);
stxt = sprintf('IQR = %5.1f/%5.1f, median = %5.1f, mean =%5.1f',...
                xlw,xup,xmd,mean(MHD));

text(0.1,0.1,stxt);
set(gca,'visible','off');

btx = 'anls_mhd_nemo_random.m';
bottom_text(btx,'pwd',1);
