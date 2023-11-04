% Plot initial MHD from SSH LCLCE contours
% to see initial errors
% Experiments with time shifted Initial Fields:  Experiments 1b
% Before analysis:
% extract LC/LCE contour
%  extr_lc_hycomPrdct_fcst.m
% compute MHD
% mhd_LCLCEcntrPrdct1b_hycom0_fcsthycom.m 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

ifcst1=10;
ifcst2=16;
irun1 = 2;
irun2 = 5;
itime1= 1;
itime2= 2;

pthmat = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
btx = 'plot_initialMHD_Prdct1b.m';

ii=0;
for iFcst=ifcst1:ifcst2
  ii=ii+1;
  for itime=itime1:itime2  % 2011 and 2012 time windows
    irr=0;
    for irun=irun1:irun2
      fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
      fmat1 = sprintf('%sMHD_LCLCE_Prdct1b_%s.mat',pthmat,fcstname);

      fprintf('Loading %s\n',fmat1);
      A=load(fmat1);  % MHD for f/cast and persistence

      
      MHD(ii).TIME(itime).MHD(irun)    = A.MHD(1,1); % initial MHD 
      MHD(ii).TIME(itime).MHD(1)       = A.MHD(1,2); % persistence
      MHD(ii).TIME(itime).iFcst(irun)  = iFcst;    % Forecast group # (groupped by hindcast ini fields)
      MHD(ii).TIME(itime).Fname        = sprintf('fcst%2.2i-%2.2i',iFcst,itime); 
    end
  end
end

FCST = [ifcst1:ifcst2];
RUN  = [irun1:irun2];
Nfgr = length(RUN);

% Bar diagram
% Choose bin width such that the bar groups  
% are separated by at least eps0
eps0=0.1; % minimum distance between the groups

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

DSHIFT = [-2; -1; 1; 2];
% Colors for different perturbed runs within same f/cast
CLR = [0   0    0;...
       0.8 0.2  0;...
       1   0.6  0; ...
       0   1    0.5; ... 
       0   0.5  1];

figure(1); clf;
set(gcf,'Position',[1428         706         874         511]);

for itime=1:2
  if itime==1
    axes('Position',[0.09 0.6 0.7 0.32]);
    TPeriod='May 2011';
  else
    axes('Position',[0.09 0.1 0.7 0.32]);
    TPeriod='Feb 2012';
  end
  hold on;
 
  ii=0;
  for iFcst=ifcst1:ifcst2
    ii=ii+1;
    for irun=irun1:irun2
      mhd = MHD(ii).TIME(itime).MHD(irun);
      clr = CLR(irun,:);

      ixx = xbin(irun-1)+ii;
      patch([ixx-dx ixx-dx ixx+dx ixx+dx],[0 mhd mhd 0],clr,'Edgecolor','none');

      Fnms{ii}=MHD(ii).TIME(itime).Fname;
    end
  end


  set(gca,'tickdir','out',...
								'xlim',[1+xbin(1)-2*dx ii+xbin(end)+2*dx],...
								'xtick',[1:ii],...
        'ytick',[0:2:20],...
        'ylim',[0 14],...
								'ygrid','on',...
								'xticklabel',Fnms,...
								'Fontsize',11);

  yr=2011;
  if itime==2, yr=2012; end
  stl = sprintf('SSH LCLCE MHD(km) Initial by f/cast groups, %i',yr);
  title(stl);
end

    
% Legend
nfcst = ifcst2-ifcst1+1;
ifc=0;
for irun=irun1:irun2  % forecast groups
		dt=DSHIFT(irun-1);
		if dt<0
				sxx = sprintf('-%iday',abs(dt));
		else
				sxx = sprintf('+%iday',abs(dt));
		end
  HName{irun}=sxx;
end
axes('Position',[0.8 0.4 0.15 0.5]);
hold on

dxx=0.2;
dyy=0.2;
ymx=(dyy+0.08)*(nfcst-1);
ifc=0;
for irun=irun1:irun2
  ifc=ifc+1;
  clr1 = CLR(irun,:);
  x1=1;
  x2=x1+dxx;
  y1=ymx-(dyy+0.08)*(ifc-1);
  y2=y1+dyy;
  patch([x1 x1 x2 x2],[y1 y2 y2 y1],clr1,'Edgecolor','none');
  stl=HName{irun};
  text(x2+dxx*0.2,(y1+y2)/2,stl,'Fontsize',12);
end
axis('equal');
set(gca,'xlim',[1 3],...
        'ylim',[0 ymx+dyy],...
        'visible','off');

bottom_text(btx,'pwd',1);





