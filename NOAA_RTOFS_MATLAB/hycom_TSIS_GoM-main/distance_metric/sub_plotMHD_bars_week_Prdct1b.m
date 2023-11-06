% Plot MHD in bars form Prdct1b experiments
% by weeks
% Group by day shift perturbations
% Show individual f/casts
%
function sub_plotMHD_bars_week_Prdct1b(nfg,CLR,POOL,Nfgr,ifcst1,ifcst2,...
                                        irun1,irun2,WK,anls_nm);

% Choose bin width such that the bar groups  
% are separated by at least eps0
eps0=0.1; % minimum distance between the groups

nruns = irun2-irun1+1;
% Half-width => dx
dx = 0.5*(1-eps0)/nruns;
% Make sure all bins < 1:
smm = nruns*dx*2+eps0;
if smm>1
  dlt = smm-1;
  dltx = dlt/(nruns*2);
  dx = dx-dltx;
end

xbin = [-0.5+eps0/2+dx:2*dx:0.5];

DSHIFT = [-2; -1; 1; 2];
nwk=length(WK)-1;

%keyboard
% Split into f/casts:
dfcst=2*dx/Nfgr;
fcst_bins=[-dx:dfcst:dx];

yl1=[];
yl2=[];
if isfield(POOL,'ylim')
  yl1=POOL(1).ylim(1);
  yl2=POOL(1).ylim(2);
else
  yl1=0;
  yl2=0;
  for ifc=1:Nfgr
    ntime = length(POOL(ifc).TIME);
    for it=1:ntime
      for irun=irun1:irun2
        dmm=POOL(ifc).TIME(it).RUN(irun).mhdPrst;
        yl2=max([yl2,max(dmm)]);
      end
    end
  end
%fprintf('yl2=%6.4g\n',yl2);
end

HName = [];

if ~isempty(nfg) 
  figure(nfg); clf;
  set(gcf,'Position',[1428         706         874         511]);
else
  figure(2); clf;
end

keyboard

for itime=1:2 % summary for 2 time period of forecasts
  if itime==1
    axes('Position',[0.07 0.6 0.8 0.32]);
    TPeriod='May 2011';
  else
    axes('Position',[0.07 0.1 0.8 0.32]);
    TPeriod='Jan 2012';
  end
  hold on;

	ntime = length(POOL(1).TIME);
	if itime>ntime, continue; end;

  for irun=irun1:irun2
		for iwk=1:nwk % 30-day intervals

      id1=WK(iwk);
      id2=WK(iwk+1)-1;
      
      ifc=0;
      mRun = [];
      DY   = [];
      DPS  = [];
      mPs  = [];
			for iFcst=ifcst1:ifcst2  % forecast groups
				ifc=ifc+1;

        mhd = POOL(ifc).TIME(itime).RUN(irun).mhd;
        id2=min([id2,length(mhd)]);
        dmm=mhd(id1:id2);

        mFc = median(dmm);
        mRun = [mRun;dmm];

        DY(ifc)=mFc;   % individual f/casts for shift irun

%  Persistence
        pmm = POOL(ifc).TIME(itime).RUN(irun).mhdPrst(id1:id2);
        DPS(ifc) = median(pmm);
        mPs = [mPs;pmm];

    % Range:
        Llw(ifc) = prctile(dmm,10);
        Lup(ifc) = prctile(dmm,90);

      end

% Overall median/mean
      mDY=median(mRun);
      clr0 = [0.9 0.9 0.9];
      ixx=xbin(irun-1)+iwk;
%
% Individual f/casts:
%      for ii=1:length(DY)
%        dy = DY(ii);
%								clr1 = CLR(ii,:);
%        plot([ixx-dx ixx+dx],[dy dy],'Color',clr1,'Linewidth',4);
%      end
% Individual f/casts:
      for ii=1:length(DY)
        dy = DY(ii);
        clr1 = CLR(ii,:);
        i1=ixx+fcst_bins(ii);
        i2=ixx+fcst_bins(ii+1);
        vx=[i1,i1,i2,i2];
        vy=[0,dy,dy,0];
        fill(vx,vy,clr1,'edgecolor','none');
% Range - percentiles:
        llw = Llw(ii);
        lup = Lup(ii);

        ix0 = 0.5*(i1+i2);
        plot([ix0 ix0],[llw lup],'-','Color',[0 0 0],'Linewidth',1.5);

      end

% Persistence:
%      for ii=1:length(DPS)
%        dy = DPS(ii);
%        clr1 = CLR(ii,:);
%        pp1=plot([ixx-0.8*dx ixx+0.8*dx],[dy dy],'--','Color',clr1,'Linewidth',2);
%      end
% Persistence:
      for ii=1:length(DPS)
        dy = DPS(ii);
        clr1 = CLR(ii,:);
%        i1=ixx+fcst_bins(ii);
%        i2=ixx+fcst_bins(ii+1);
        i1=ixx-dx;
        i2=ixx+dx;
        plot([i1 i2],[dy dy],'--','Color',clr1,'Linewidth',1.5);
      end
        


% Label of perturb day shifts:
						ixx=xbin(irun-1)+iwk-dx;
						dt=DSHIFT(irun-1);
						if dt<0
								sxx = sprintf('-%i',abs(dt));
						else
								sxx = sprintf('+%i',abs(dt));
						end
						text(ixx,-5,sxx,'Fontsize',8);
%keyboard


    end
  end
  set(gca,'tickdir','out',...
          'xlim',[1+xbin(1)-2*dx iwk+xbin(end)+2*dx],...
          'xtick',[1:15],...
          'ygrid','on',...
          'Fontsize',12);
%  dstr1=MHD(j1).Date_str;
%  dstr2=MHD(j2).Date_str;

  if ~isempty(yl1)
    set(gca,'ylim',[yl1 yl2]);
    if yl2<=100
      set(gca,'ytick',[0:10:200]);
    elseif yl2>100 & yl2<=200
      set(gca,'ytick',[0:10:200]);
    elseif yl2>200
      set(gca,'ytick',[0:50:200]);
    elseif yl2<=0.1
      set(gca,'ytick',[0:0.01:1]);
    end
  end

% Separate weeks:
  for iwk=1:nwk-1
    plot([iwk+0.5 iwk+0.5],[0 yl2],'--','color',[0.7 0.7 0.7]);
  end;


  stl=sprintf('%s Init:%s',anls_nm,TPeriod);
  title(stl);
end


nfcst = ifcst2-ifcst1+1;
ifc=0;
for iFcst=ifcst1:ifcst2  % forecast groups
  ifc=ifc+1;
  HName{ifc}=sprintf('F/cast %2.2i',iFcst);
end
axes('Position',[0.88 0.4 0.11 0.5]);
hold on

dxx=0.2;
dyy=0.2;
ymx=(dyy+0.08)*(nfcst-1);
for ifc=1:nfcst
  clr1 = CLR(ifc,:);
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
