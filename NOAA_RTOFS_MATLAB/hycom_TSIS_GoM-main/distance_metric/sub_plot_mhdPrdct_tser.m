function sub_plot_mhdPrdct_tser(MHD,CLR,irun,ifn,ipst,itot,ifcst);
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
clr = CLR(irun,:);
mhd = MHD(ifcst).mhd;
plot(mhd,'-','Color',clr,'Linewidth',2);
ymx=max([max(mhd,ymx)]);

clr=CLR(ipst,:);
mhd=MHD(ifcst).nemo0;
plot(mhd,'-','Color',clr,'Linewidth',2);
ymx=max([max(mhd,ymx)]);

ylm2=1.1*ymx;
nrc = length(mhd);
set(gca,'tickdir','out',...
								'xlim',[1 nrc],...
								'ylim',[0 ylm2],...
								'xtick',[0:10:100],...
								'xgrid','on',...
								'ygrid','on',...
								'Fontsize',12);
%keyboard
dstr = MHD(ifcst).Date_str;
%if itime==1
%		dstr = sprintf('May/June 2011');
%else
%		dstr = sprintf('Jan/Feb 2012');
%end;
%itot = (nhnd-7+itime-1)*Nruns+irun;
stl = sprintf('MHD(km) HYCOM %s vs NEMO, LCLCE cntrs, %s ',MHD(ifcst).Name,dstr);
title(stl,'Fontsize',12);

xlabel('Fcst days');
ylabel('MHD, km');


%
%  Cumulative score
%
ymx=0;
axes('Position',[0.08 0.1 0.85 0.32]);
hold on;
clr = CLR(irun,:);
mhd = MHD(ifcst).mhd;
plot(cumsum(mhd),'-','Color',clr,'Linewidth',2);
ymx=max([max(sum(mhd),ymx)]);

clr=CLR(ipst,:);
mhd=MHD(ifcst).nemo0;
plot(cumsum(mhd),'-','Color',clr,'Linewidth',2);
ymx=max([max(sum(mhd),ymx)]);

ylm2=1.1*ymx;
set(gca,'tickdir','out',...
								'xlim',[1 nrc],...
								'ylim',[0 ylm2],...
								'xtick',[0:10:100],...
								'xgrid','on',...
								'ygrid','on');
title('Cum score');


% Legend:
axes('Position',[0.08 0.7 0.4 0.25]);
hold;
x0=0;
clr = CLR(irun,:);
nm = MHD(irun).Name;
plot([0 0.5],[2 2],'-','Color',clr,'Linewidth',2);
text(0.6,2,nm,'Fontsize',12);
clr = CLR(ipst,:);
nm = 'Persistence';
plot([0 0.5],[1 1],'-','Color',clr,'Linewidth',2);
text(0.6,1,nm,'Fontsize',12);


set(gca,'xlim',[-0.2 4], ...
								'ylim',[-2 4],...
								'xtick',[],...
								'ytick',[],...
								'visible','off');



return
