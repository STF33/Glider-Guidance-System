% Plot SST from NEMO
% 30th and 60th NEMO grid points
% and UGOS PIES locations
% for OSSEs
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

%pthdat = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/obs/aviso/'; % for 2009
%pthdat  = '/nexsan/people/abozec/TSIS/data/aviso/';  % <-- all data here, gzip
pthdat = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/data/aviso/';  % copied from Alex's dir
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';

% In TSIS - combine 4 days of altimetry to get more data
%dnmb = datenum(2011,06,01);
dnmb0  = datenum(2011,6,1);  % analysis date:
dnmb1  = dnmb0-4;
dnmb2  = dnmb0+4;

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

ln1=min(min(LON));
ln2=max(max(LON));
lt1=min(min(LAT));
lt2=max(max(LAT));

% Colormap for SST field:
% Colormap
ncc=50;
cl1=[0.8,0.2,0.8];
cl2=[0.6,0.2,0.9];
clr1=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[0.,0.4,0.9];
clr2=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[0,0.9,1];
clr3=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[0,0.7,0.4];
clr4=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[1,1,0.3];
clr5=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[0.9,0.4,0];
clr6=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[.9,0.1,0];
clr7=mix_2colors(cl1,cl2,ncc);

cmp_sst = [clr1;clr2;clr3;clr4;clr5;clr6;clr7];
cmp_sst = smooth_colormap(cmp_sst,5);

pthnemo = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
fprintf('Loading NEMO grid %s\n',fgrd);
load(fgrd);

% Det 30th and 60th NEMO indices
% Deep ocean > 20 m
z0 = -25;
iz0 = max(find(ZZN>=z0));
Tz0 = sub_get_NEMO_TS(dnmb0,fldnm,iz0);
% Only GoM:
Tz0(LONN>-81) = NaN;
Tz0(LONN>-90 & LATN<21) = NaN;
Tz0(LONN>-85 & LATN<22) = NaN;
[mN,nN] = size(Tz0);
[II,JJ]=meshgrid([1:nN],[1:mN]);
di=30;
sX  = LONN(1:di:end,1:di:end);
sY  = LATN(1:di:end,1:di:end);
sT  = Tz0(1:di:end,1:di:end);
I30 = find(~isnan(sT)); 
sX30 = sX(I30);
sY30 = sY(I30);

di=60;
sX  = LONN(1:di:end,1:di:end);
sY  = LATN(1:di:end,1:di:end);
sT  = Tz0(1:di:end,1:di:end);
I60 = find(~isnan(sT)); 
sX60 = sX(I60);
sY60 = sY(I60);


% Read SST:
fldnm = 'toce';
iz0=1;
TT = sub_get_NEMO_TS(dnmb0,fldnm,iz0);

figure(5); clf;
set(gcf,'Position',[1491         726         791         583]);
axes('Position',[0.08 0.2 0.85 0.72]);
hold on;

pcolor(LONN,LATN,TT); shading flat;
colormap(cmp_sst);
ct1=25;
ct2=29;
caxis([ct1 ct2]);

% Every 30th NEMO point:
f_di=60;
if f_di==30
  plot(sX30,sY30,'k.');
elseif f_di==60
  plot(sX60,sY60,'k.');
end

sttl = sprintf('NEMO SST, dx=%i, analysis date=%s',f_di,datestr(dnmb0));  

% Add UGOS PIES (=1) or UGOS extended PIES (=2)
pthosse = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/OSSE/';
YR=2011;
f_pies = 0;
if f_pies==1
  finp=sprintf('%sugos_mooring_%i.mat',pthosse,YR);
  fprintf('Loading %s\n',finp);
  A=load(finp);
  plat=A.ugos_lat;
  plon=A.ugos_lon;
  sttl = sprintf('NEMO SST, PIES, analysis date=%s',...
     datestr(dnmb0));
elseif f_pies==2
  finp=sprintf('%sextd_mooring_2011.mat',pthosse);
  A=load(finp);
  plat=A.extd_lat;
  plon=A.extd_lon;
  sttl = sprintf('NEMO SST, extd PIES, analysis date=%s',...
     datestr(dnmb0));
end
if f_pies>0
  plot(plon,plat,'^','Color',[0 0 0],'linewidth',1.5);
end

title(sttl);

axis('equal');
set(gca,'tickdir','out',...
    'xlim',[ln1 -80],...
    'ylim',[18 31.2]);
set(gca,'Color',[0 0 0]);

%
% Plot colorbar
%
clb=colorbar('SouthOutside');

set(clb,'Position',[0.1 0.1 0.66 0.025],...
        'Fontsize',13,...
        'Ticks',[10:1:40],...
        'Ticklength',0.025);

btx = 'plot_NEMOsst_PIES.m';
bottom_text(btx,'pwd',1);




