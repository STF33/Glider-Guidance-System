% Plot MHD in bars with 10-90 prctiles
% by weeks - sort from best to worst fcst
%
function sub_plotMHD_barsrt_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,WK,anls_nm);

nwk=length(WK)-1;  % # of forecast weeks

% Choose bin width such that the bar groups  
% are separated by at least eps0
eps0=0.2;

% Half-width => dx
dx = 0.5*(1-eps0)/Nfgr;
% Make sure all bins < 1:
smm = Nfgr*dx*2+eps;
if smm>1
  dlt = smm-1;
  dltx = dlt/(Nfgr*2);
  dx = dx-dltx;
end

xbin = [-0.5+eps0/2+dx:2*dx:0.5];


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
      dmm=POOL(ifc).Time(it).pm3;
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
    TPeriod='Apr/May 2011';
  else
    axes('Position',[0.09 0.1 0.7 0.32]);
    TPeriod='Jan/Feb 2012';
  end
  hold on;

  for ifc=1:Nfgr  % forecast groups
    iFcst=IFCST(ifc);
    clr1 = CLR(iFcst,:);

    ntime = length(POOL(ifc).Time);
    if itime>ntime, continue; end;

    for iwk=1:nwk %
%keyboard
      dmm=POOL(ifc).Time(itime).WK(iwk).pw;
      if isempty(dmm); continue; end;

      mFc = median(dmm);
      Fc1 = prctile(dmm,25); % interquartiles
      Fc2 = prctile(dmm,75);
      SC(ifc,iwk) = mFc;
      FC1(ifc,iwk) = Fc1;
      FC2(ifc,iwk) = Fc2;
    end
  end

  for iwk=1:nwk
% Rank by cum score:
    [BB,I] = sort(SC(:,iwk),'ascend');
   
    for ix0=1:length(I)
      ifc=I(ix0);
      iFcst=IFCST(ifc);
      clr1 = CLR(iFcst,:);

      ixx=xbin(ix0)+iwk;
      dy=SC(ifc,iwk);
      patch([ixx-dx ixx-dx ixx+dx ixx+dx],[0 dy dy 0],clr1,'Edgecolor','none');
    % Range:
      llw = FC1(ifc,iwk);
      lup = FC2(ifc,iwk);
      plot([ixx-0.03 ixx+0.03],[llw llw],'k-');
      plot([ixx-0.03 ixx+0.03],[lup lup],'k-');
      plot([ixx ixx],[llw lup],'k-');
%keyboard
%      if iwk==13 & ifc==6; keyboard; end;
      if Fc1>mFc | Fc2<mFc,
        fprintf('*** ERR: Check median and IQrange \n');
        keyboard;
      end
    end
  end
%
% Separate weeks:
  for iwk=1:nwk-1
    plot([iwk+0.5 iwk+0.5],[0 yl2],'--','color',[0.7 0.7 0.7]);
  end;

  set(gca,'tickdir','out',...
          'xlim',[0.5 nwk+0.6],...
          'xtick',[1:nwk],...
          'ygrid','on',...
          'Fontsize',14);
  if ~isempty(yl1)
    set(gca,'ylim',[yl1 yl2]);
    if yl2<=100 & yl2>40.
      set(gca,'ytick',[0:10:200]);
    elseif yl2>100 & yl2<=200
      set(gca,'ytick',[0:20:200]);
    elseif yl2>200
      set(gca,'ytick',[0:50:200]);
    elseif yl2<=0.1
      set(gca,'ytick',[0:0.01:1]);
    elseif yl2>0.02 & yl2<0.2
      set(gca,'ytick',[0:0.025:1]);
    elseif yl2>0.05 & yl2<0.5
      set(gca,'ytick',[0:0.05:1]);
    end
  end

  stl=sprintf('%s Init:%s',anls_nm,TPeriod);
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

%  if iHnd>=10, iHnd=ifc; end;   % for 1st set of f/casts (<10) use colors matching hindcasts

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
