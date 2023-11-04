% Plot all RMSE
% time series
%
function sub_plotRMSEall_tser(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);

yl1=[];
yl2=[];
if isfield(POOL,'ylim')
  yl1=POOL(1).ylim(1);
  yl2=POOL(1).ylim(2);
else
  yl1=0;
  yl2=0;
  for ifc=1:Nfgr
    ntime = length(POOL(ifc).Time);
    for it=1:ntime
      dmm=POOL(ifc).Time(it).tser;
      yl2=max([yl2,max(dmm)]);
    end
  end
end

HName = [];

figure(nfg); clf;
set(gcf,'Position',[1428         706         874         511]);

for itime=1:2 % summary for 2 time period of forecasts
  if itime==1
    axes('Position',[0.09 0.6 0.7 0.32]);
    TPeriod='May 2011';
  else
    axes('Position',[0.09 0.1 0.7 0.32]);
    TPeriod='Jan 2012';
  end
  hold on;

  for ifc=1:Nfgr  % forecast groups
    iFcst=IFCST(ifc);
% For 1st set of experiments - match colors with HYCOM hindcasts
% used for initialization
    if (iFcst<10) 
      clr1 = CLR(iFcst,:);
    else
      clr1 = CLR(ifc,:);
    end

    ntime = length(POOL(ifc).Time);
    if itime>ntime, continue; end;

    rmse = POOL(ifc).Time(itime).tser;
    tt = [0:length(rmse)-1]; 
    plot(tt,rmse,'Linewidth',2,'Color',clr1);
  end

  if ~isempty(yl1)
    set(gca,'ylim',[yl1 yl2]);
    if yl2>100
      set(gca,'ytick',[0:20:200]);
    elseif yl2<=0.1
      set(gca,'ytick',[0:0.01:1]);
    end
  end

  set(gca,'tickdir','out',...
          'xlim',[0 length(tt)],...
          'xtick',[0:10:120],...
          'xgrid','on',...
          'ygrid','on',...
          'Fontsize',14);

  stl=sprintf('%s Init: %s',anls_nm,TPeriod);
  title(stl);
end
for ifc=1:Nfgr  % forecast groups
  iFcst=IFCST(ifc);
  HName{ifc}=EXPT(iFcst).Name_short;
end
axes('Position',[0.8 0.5 0.15 0.4]);
hold on


dxx=0.2;
dyy=0.2;
ymx=(dyy+0.08)*(Nfgr-1);
for ifc=1:Nfgr
  iHnd = IFCST(ifc);

  if iHnd>=10, iHnd=ifc; end;   % for 1st set of f/casts (<10) use colors matching hindcasts

  clr1 = CLR(iHnd,:);
  x1=1;
  x2=x1+dxx;
  y1=ymx-(dyy+0.08)*(ifc-1);
  y2=y1+dyy;
  patch([x1 x1 x2 x2],[y1 y2 y2 y1],clr1,'Edgecolor','none');
  stl=HName{ifc};
  text(x2+dxx*0.2,(y1+y2)/2,stl,'Fontsize',12);
end
axis('equal');
set(gca,'xlim',[1 3],...
        'ylim',[0 ymx+dyy],...
        'visible','off');

return
