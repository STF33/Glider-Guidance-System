% Plot LC/LCE and ssh
% OSE hindcasts 
% For specified date

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

yr1 = 2009;
yr2 = 2009;
%esim = 'noPIES';
%esim = 'PIES';
esim = 'PIESv2';   % updated hindcast 
Z0 = -200;
f_save = 318;    % make f_save > 0 to save fig and skip N<f_save fields

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';
pthfig  = sprintf('/Net/mars/ddmitry/hycom/hycom_TSIS/frames_ssh_%s/',esim);

% HYCOM:
rg=9806;  % convert pressure to depth, m
huge=1e20;
Bisol = 0.17;  % ssh contour

%Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
HH(isnan(HH))=100;


% GoM region HYCOM:
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
INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
Iocn = find(INH==1 & HH<Z0);
Nocn=length(Iocn);
clear XM YM

Imean = INH;
Imean(HH>Z0)=0;
% 
% Subsample to a smaller domain:
xnas1 = min(LON(Iocn));
xnas2 = max(LON(Iocn));
ynas1 = min(LAT(Iocn));
ynas2 = max(LAT(Iocn));
[J,I]=find(LON<xnas1);
its1=max(I);
[J,I]=find(LON>xnas2);
its2=min(I);
[J,I]=find(LAT<ynas1);
jts1=max(J);
[J,I]=find(LAT>ynas2);
jts2=min(J);

xt1=LON(jts1,its1);
yt1=LAT(jts1,its1);
xt2=LON(jts2,its2);
yt2=LAT(jts2,its2);


fmatout = sprintf('%shycom_LCcontour_OSEhindcast_%s_2009-2011.mat',...
pthmat,esim);

if strncmp(esim,'PIESv2',6)
  fmatout=sprintf('%shycom_LCcontour_OSEhindcast_%s_2009-2009.mat',...
pthmat,esim);
end

dmean=1;




YRPLT=[];
cc=0;
for iyr=yr1:yr2
  id1=1;
  id2=365;
  ic=mod(iyr,4);
  if ic==0 & id2==365, id2=366; end;
  if iyr>yr1
    id1=1;
  end
  for iday=id1:id2
    cc=cc+1;
    jd1=datenum(iyr,1,1);
    dnmb=jd1+iday-1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=iday;
    YRPLT(cc,3)=dnmb;
    DV(cc,:) = datevec(dnmb);
  end
end
nrc=size(YRPLT,1);


% Colormap
ncc=100;
cl1=[0,0.3,0.6];
cl2=[1,1,1];
clr1=mix_2colors(cl1,cl2,ncc);

%cl1=cl2;
%cl2=[0.6,1,0.8];
cl1=cl2;
cl2=[0.6,0.9,0.7];
clr2=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[1,0.8,0.2];
clr3=mix_2colors(cl1,cl2,ncc);
cmp = [clr1;clr2;clr3];
cmp = smooth_colormap(cmp,5);


ifg=0;
for ii=1:nrc
  tic;
  yr  = DV(ii,1);
  mo  = DV(ii,2);
  dm  = DV(ii,3);
  dnmb= YRPLT(ii,3);
  iday= YRPLT(ii,2);
  HR  = 0;

  switch (esim)
   case({'PIES','noPIES'})
    pthi=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_%s/',yr,esim);
   case('PIESv2')
    pthi='/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_newtsis/gofs30_withpies/';  
  end
  fina = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthi,yr,iday,HR);
  finb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthi,yr,iday,HR);

  if ~exist(fina,'file');
    fprintf('Not found %s\n',fina);
    continue;
  end

  ifg = ifg+1;
  if ii<f_save; 
    fprintf('Skipping day %s\n',datestr(dnmb));
    continue;
  end

  ssh = sub_getSSH_hycom(fina,finb,Imean,dmean);
  sshD=ssh;
  sshD(INH==0)=nan;

%
%  LC/LCE contours:
  fprintf('Loading %s\n',fmatout);
  load(fmatout); 
% Find time index:
  TMLC = LCXY.TM;
  TMLCE = LCE(1).TM;
  ilc = find(TMLC==dnmb);
  ilce = find(TMLCE==dnmb);

  btx = 'frames_ssh_LCcntr_OSEhindcast.m';
  figure(1); clf;
  set(gcf,'Position',[1737         707         760         604]);
  pcolor(LON,LAT,sshD); shading flat;
  hold on;
  contour(LON,LAT,sshD,[Bisol Bisol],'k');
  contour(LON,LAT,HH,[0 0],'Color',[0.5,0.5,0.5]);
  colormap(cmp);

  clm1=-0.4;
  clm2=0.8;
  dcc=0.2;
  caxis([clm1 clm2]);

  axis('equal');
  set(gca,'xlim',[-98 -80],...
      'ylim',[17 31]);

  xlc = LCXY.XY(ilc).X;
  ylc = LCXY.XY(ilc).Y;
  dstr = datestr(dnmb);
  plot(xlc,ylc,'r.');

  nlce = length(LCE);
  for jj=1:nlce
    nxy = length(LCE(jj).XY);
    if ilce>nxy, continue; end;
    xlce = LCE(jj).XY(ilce).X;
    ylce = LCE(jj).XY(ilce).Y;
    if isempty(xlce), continue; end;

    plot(xlce,ylce,'r.');
  end
  stl = sprintf('%s %s',esim,dstr);
  title(stl);
  ax1 = gca;

  axes('Position',[0.44, 0.15, 0.45, 0.08]);
%  set(gca,'visible','off');
  set(gca,'xtick',[],...
          'ytick',[],...
          'box','on');
%  ax2=gca;
%  gca=ax1;
%  hb = colorbar('peer',ax2,'south');
  hb = colorbar('south');
  caxis([clm1, clm2]);
  set(hb,'Ticks',[clm1:dcc:clm2],...
         'Fontsize',10);


  bottom_text(btx,'pwd',1);

  if f_save > 0
    fgnm = sprintf('%ssshLC_OSE_%s_%4.4i.png',pthfig,esim,ifg);
    fprintf('Saving %s\n',fgnm);
%keyboard
    print('-dpng','-r200',fgnm);
    fprintf(' Process time %7.5f min\n\n',toc/60);
  end
  
end


