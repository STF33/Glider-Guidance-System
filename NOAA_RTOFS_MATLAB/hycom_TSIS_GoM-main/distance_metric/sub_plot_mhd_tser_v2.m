function sub_plot_mhd_tser_v2(mhd,mhdP,clr1,clr2,ifn,stl,ymax,irun);
% Time series of MHD scores 
% for the forecast/hindcast experiments
figure(ifn); clf;
set(gcf,'Position',[1134 750 958 566]);

%
%   MHD time series
%

axes('Position',[0.08 0.6 0.85 0.32]);
hold on;
ymx=0;
plot(mhd,'-','Color',clr1,'Linewidth',2);
ymx=max([max(mhd,ymx)]);

plot(mhdP,'-','Color',clr2,'Linewidth',2);
ymx=max([max(mhdP,ymx)]);

nrc = length(mhd);
if ~isempty(ymax), 
  ylm2=1.1*ymax;
else
  ylm2=1.1*ymx;
end

set(gca,'tickdir','out',...
								'xlim',[1 nrc],...
								'ylim',[0 ylm2],...
								'xtick',[1:10:100],...
								'xgrid','on',...
								'ygrid','on',...
								'Fontsize',12);

%stl = sprintf('MHD(km) HYCOM Frcst %s vs NEMO, LCLCE cntrs, %s ',MHD(itot).Name,dstr);
title(stl,'Fontsize',12);

xlabel('Fcst days');
ylabel('MHD, km');

% No legend
if irun<=0, return; end; 

% Legend:
axes('Position',[0.1 0.2 0.3 0.3]);
hold on;
x1=0.1;
x2=x1+0.2;
y1=1;

y1=y1-0.1;
clr=clr1;
plot([x1 x2],[y1 y1],'-','Color',clr,'Linewidth',2);
stxt = sprintf('run=%i ',irun);
text(x2+0.1,y1,stxt,'Fontsize',12);

y1=y1-0.1;
clr=clr2;
plot([x1 x2],[y1 y1],'-','Color',clr,'Linewidth',2);
stxt = 'Persist';
text(x2+0.1,y1,stxt,'Fontsize',12);
set(gca,'xlim',[x1 x1+0.8],...
								'ylim',[y1-0.1 1], ...
								'Visible','off');

return
