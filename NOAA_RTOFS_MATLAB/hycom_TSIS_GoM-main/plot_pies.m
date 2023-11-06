% Plot ssh subtacted area average SSH from TSIS analysis
% submit:  matlab -nodesktop -nosplash < plot_ssh5.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

% ----------------
% Flags
% ---------------
s_fig = 1;

ys=2009;
ye=2010;
dday = 1; % day stepping for plotting
ds=1;
%ms=305;
%dJ1 = datenum(ys,1,1);
%nday0=dJ1+ds-1;

de=2;
id1=2;
id2=365;

rg=9806;  % convert pressure to depth, m
huge=1e20;

pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/frames_deepUpies/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthdat = '/Net/gleam/morey/NASGRP/';
fin = sprintf('%sdata/DynLoopData4Steve.mat',pthdat);
load(fin);
PS=dynloop;
TM=datenum(2009,1,1,0,0,0)+PS.dd;
npp = length(TM);

plat = PS.lat;
plon = PS.lon;

%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/GoM/GOMl0.04/%s/data_anls/',expt);
btx = 'plot_pies.m';

% Read the topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);

Lmsk = HH*0;
Lmsk(HH<0) = 1;


% GoM region:
GOM=[366   489
   476   531
   583   560
   576   646
   508   827
   336   848
   204   829
    64   798
    19   746
    16   662
    12   578
    25   455
    71   382
   165   356
   281   400];

[XM,YM]=meshgrid([1:n],[1:m]);
IN = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
clear XM YM

IDp = find(IN==1 & HH<-1000);
IB  = find(IN==0 | HH>-1000);

HH(IB)=nan;

cmpL = [0 0 0; .4 .4 .4];

c1=0.;
c2=0.3;
CMP = create_colormap_WBYR(200,c1,c2);
cmp = CMP.colormap;
cnt = CMP.intervals;
nint= length(cnt);
%
cntr = 0;
ff=figure('Visible','off');
for it=2:2:npp
  tic;
  dnmb = TM(it);
  dv=datevec(dnmb);
  fprintf('Plotting %4.4i/%2.2i/%2.2i %2.2i:%2.2i\n',dv(1:5));

  Ub = squeeze(PS.u3000(:,:,it))*0.01; % cm/s->m/s
  Vb = squeeze(PS.v3000(:,:,it))*0.01; % cm/s->m/s
  Sb = sqrt(Ub.^2+Vb.^2);
% ==========================
%   Plot deep flow:
% ==========================
  clf;
%      contour(lon,lat,HH,[-100 -100],'Color',[0.7 0.7 0.7],'linewidth',1);
  pcolor(LON,LAT,Lmsk); shading flat;
  colormap(cmpL);
  freezeColors;

  hold on;

  pcolor(plon,plat,Sb); shading flat;
  caxis([c1 c2]);
  colormap(cmp);

  contour(LON,LAT,HH,[-100 -100],'Color',[0.8 0.8 0.8]);
  contour(LON,LAT,HH,[-5000:1000:-1000],'Color',[0.6 0.6 0.6]);

  axis('equal');
  set(gca,'tickdir','out',...
	  'xlim',[-97.5 -83.8],...
	  'ylim',[19 30]);

% Plot vectors:
  fprintf('Drawing vectors ...\n');
  Ib=find(~isnan(Ub));
  dltx=5;
  cf=0.25;
  beta=30;
  v_col=[0 0. 0];
%    v_col=[1 0 0];
  lwd=1.;
  scl=0.25;
  for iib=1:dltx:length(Ib)
    I0=Ib(iib);
    [j,i] = ind2sub(size(Ub),I0);
    if ~isnan(Sb(j,i)) & Sb(j,i)>0.02
      u=Ub(j,i);  % m/s
      v=Vb(j,i);
      ss=sqrt(u.^2+v.^2);
      u=u/ss;
      v=v/ss;
      x1=plon(j,i);
      x2=x1+u*scl;
      y1=plat(j,i);
      y2=y1+v*scl;
      draw_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
    end;
  end;  % indices loop

% Plot ssh
%  contour(LON,LAT,ssh,[0:0.1:0.9],'Color',[1 1 0.7]);
%  contour(LON,LAT,ssh,[-0.5:0.05:-0.01],'k--','Color',[1 1 0.7]);

%  title(['SSH, ',regn,' - ',expt,' Year: ',int2str(year),'; D= ',sday]);
  date_str = sprintf('%4.4i/%2.2i/%2.2i %2.2i:%2.2i',dv(1:5));
  stt=sprintf('PIES U, 3000m, %s',date_str);
  title(stt,'Fontsize',11);
  

  chb = colorbar;
  set(chb,'Position',[0.91 0.12 0.02 0.78],...
	  'TickLength',0.032);

%  bottom_text(btx,'pwd',1,'Fontsize',4);
%  drawnow

  cntr=cntr+1;  
  if s_fig==1
%    fnm=sprintf('goml%s-SSH%4.4i_%s',expt,year,sday);
    fnm=sprintf('pies_deepU-%4.4i',cntr);
    fout=[pthfig,fnm];
    fprintf('Saving %s\n',fout);
    set(gcf,'InvertHardcopy','off','Color',[1 1 1]);
    print('-dpng','-r250',fout);
  end
  
  fprintf('1 step processed: %6.3f min\n\n',toc/60);
%keyboard
end;


