% Plot LC/LCE and ssh
% OSE hindcasts 
% For specified date

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

dnmb = datenum(2009,09,02);
%esim = 'noPIES';
%esim = 'PIES';
esim = 'PIESv2';   % updated hindcast 
Z0 = -200;

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';


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
DV = datevec(dnmb);
yr  = DV(1);
mo  = DV(2);
dm  = DV(3);
iday= dnmb-datenum(yr,1,1)+1; 
HR  = 0;

switch (esim)
 case({'PIES','noPIES'})
  pthi=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_%s/',yr,esim);
 case('PIESv2')
  pthi='/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_newtsis/gofs30_withpies/';  
end

fina = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthi,yr,iday,HR);
finb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthi,yr,iday,HR);

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


btx = 'plot_ssh_LCcntr_OSEhindcast.m';
figure(1); clf;
pcolor(LON,LAT,sshD); shading flat;
hold on;
contour(LON,LAT,sshD,[Bisol Bisol],'k');

axis('equal');
set(gca,'xlim',[-98 -80],...
				'ylim',[17 31]);

xlc = LCXY.XY(ilc).X;
ylc = LCXY.XY(ilc).Y;
dstr = datestr(dnmb);
plot(xlc,ylc,'r.');

nlce = length(LCE);
for jj=1:nlce
 xlce = LCE(jj).XY(ilce).X;
 ylce = LCE(jj).XY(ilce).Y;
 if isempty(xlce), continue; end;

 plot(xlce,ylce,'r.');
end
stl = sprintf('%s %s',esim,dstr);
title(stl);
bottom_text(btx,'pwd',1);


