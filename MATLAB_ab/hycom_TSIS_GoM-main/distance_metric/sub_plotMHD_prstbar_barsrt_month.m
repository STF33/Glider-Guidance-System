% Plot MHD in bars with IQR
% and persistence bar - next to each experiment
% for 1-3 months - sort from best to worst f/cast for each month
% 
function sub_plotMHD_barsrt_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,PRST,anls_nm);

% Choose bin width such that the bar groups  
% are separated by at least eps0
eps0=0.15; % minimum distance between the groups

% Half-width => dx
dx = 0.5*(1-eps0)/Nfgr;
% Make sure all bins < 1:
smm = Nfgr*dx*2+eps0;
if smm>1
  dlt = smm-1;
  dltx = dlt/(Nfgr*2);
  dx = dx-dltx;
end

xbin = [-0.5+eps0/2+dx:2*dx:0.5];

%xbin=[1:Nfgr];
%xbin=(xbin-mean(xbin))*2*dx;
%
%keyboard

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
      dmm=prctile(POOL(ifc).Time(it).pm3,90);
      yl2=max([yl2,max(dmm)]);
    end
  end
%fprintf('yl2=%6.4g\n',yl2);
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
%keyboard
  clear FC1 FC2 SC
  for ifc=1:Nfgr  % forecast groups
    iFcst=IFCST(ifc);
    clr1 = CLR(iFcst,:);

    ntime = length(POOL(ifc).Time);
    if itime>ntime, continue; end;

    for jm=1:3 % 30-day intervals
      if jm==1
        dmm = POOL(ifc).Time(itime).pm1;
      elseif jm==2
        dmm = POOL(ifc).Time(itime).pm2;
      else
        dmm = POOL(ifc).Time(itime).pm3;
      end

      mFc = median(dmm);
      Fc1 = prctile(dmm,25); % interquartiles
      Fc2 = prctile(dmm,75);
      SC(ifc,jm)  = mFc;
      FC1(ifc,jm) = Fc1;
      FC2(ifc,jm) = Fc2;
    end
  end

  dbr = 0.02;
  dxx = 0.8*dx;
  for jm=1:3
    [BB,I] = sort(SC(:,jm),'ascend');

    for ix0=1:length(I)
      ifc=I(ix0);
      iFcst=IFCST(ifc);
      clr1 = CLR(iFcst,:);

      ixx=xbin(ix0)+jm;
      dy=SC(ifc,jm);
      ixM = ixx-0.5*dxx;
      patch([ixx-dxx ixx-dxx ixx ixx],[0 dy dy 0],clr1,'Edgecolor','none');
    % Range:
      llw = FC1(ifc,jm);
      lup = FC2(ifc,jm);
      plot([ixM-dbr ixM+dbr],[llw llw],'k-');
      plot([ixM-dbr ixM+dbr],[lup lup],'k-');
      plot([ixM ixM],[llw lup],'k-');

    end
  end
% Persistence:
  clrP = [0.7 0.7 0.7];
  for jm=1:3
    [BB,I] = sort(SC(:,jm),'ascend');

    for ix0=1:length(I)
      ifc=I(ix0);
      iFcst=IFCST(ifc);
      clr1 = CLR(iFcst,:);

      if jm==1
        dmm = PRST(ifc).Time(itime).pm1;
      elseif jm==2
        dmm = PRST(ifc).Time(itime).pm2;
      else
        dmm = PRST(ifc).Time(itime).pm3;
      end

      mFc = median(dmm);
      Fc1 = prctile(dmm,25); % interquartiles
      Fc2 = prctile(dmm,75);

      ixx=xbin(ix0)+jm;
      dy = mFc;
      ixM = 0.5*(ixx+(ixx+dxx));
      patch([ixx ixx ixx+dxx ixx+dxx],[0 dy dy 0],clrP,'Edgecolor','None');
%      plot([ixx ixx],[0 dy],'linewidth',2.5,'color',clr1);
%      plot([ixx+dxx ixx+dxx],[0 dy],'linewidth',2.5,'color',clr1);
%      plot([ixx ixx+dxx],[dy dy],'linewidth',2.5,'color',clr1);

    % Range:
      llw = Fc1;
      lup = Fc2;
      plot([ixM-dbr ixM+dbr],[llw llw],'k-');
      plot([ixM-dbr ixM+dbr],[lup lup],'k-');
      plot([ixM ixM],[llw lup],'k-');

    end
  end



  set(gca,'tickdir','out',...
          'xlim',[0.5 3.5],...
          'xtick',[1:3],...
          'ygrid','on',...
          'xticklabel',{'0-30','31-60','61-91'},...
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
HName{ifc+1}='Persist';
axes('Position',[0.8 0.4 0.15 0.5]);
hold on

dxx=0.2;
dyy=0.2;
ymx=(dyy+0.08)*Nfgr;
for ifc=1:Nfgr
  iHnd = IFCST(ifc);

  clr1 = CLR(iHnd,:);
  x1=1;
  x2=x1+dxx;
  y1=ymx-(dyy+0.08)*(ifc-1);
  y2=y1+dyy;
  patch([x1 x1 x2 x2],[y1 y2 y2 y1],clr1,'Edgecolor','none');
  stl=HName{ifc};
  text(x2+dxx*0.2,(y1+y2)/2,stl,'Fontsize',12);
end
y1 = ymx-(dyy+0.08)*ifc;
y2 = y1+dyy;
patch([x1 x1 x2 x2],[y1 y2 y2 y1],clrP,'Edgecolor','none');
%plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],'-','Linewidth',2,'Color',clrP);
stl=HName{ifc+1};
text(x2+dxx*0.2,(y1+y2)/2,stl,'Fontsize',12);
%keyboard

axis('equal');
set(gca,'xlim',[1 3],...
        'ylim',[0 ymx+dyy],...
        'visible','off');



return
