% Plot HYCOM interpolated T onto z level
% interpolation performed in: interp_hycom2z.m
% and extracted LC / LCE contours: extr_lc_temp_hycom.m
% 
% and NEMO fields
%
%
% For LC contour finding: interpolate T or S fields to 
% Z depth Z0
%
% Extract LC and LCEs for calculating MHD
% Using T fields (demeaned) at depth Z0, contour T anomaly
%
% NEMO data
%
% NEMO data:
%ncdump -h https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc
%fin='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc';
%ncid = netcdf.open(fin);
%vid  = netcdf.inqVarID(ncid,'ssh');
%ssh  = netcdf.getVar(ncid,vid,'single');
%[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,vid)
%varname = netcdf.inqVar(ncid,vid)
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_TSIS/interp2grid
startup;

close all
clear

Z0 = -200;  % depth
T0 = 2.5;   % T contour
irec = 173;  % day

%EXON = zeros(9,1);
%EXON(9) = 1; % select expt to plot
% Select hindcast experiment #: 2 -- 9
ixx = 8;  % hindcast/free run # expt



pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

% HYCOM-TSIS hindcast experiments:
fhnd = 'hycom_tsis_expts.mat';
load(fhnd);
Nruns = length(EXPT);



clr1 = colormap_pink(200);
for ik=1:15
  clr1(ik,:)=[1 1 1];
end
clr1 = smooth_colormap(clr1,5);
clr1 = flipud(clr1);

clr2 = colormap_turt(200);
for ik=1:15
  clr2(ik,:)=[1 1 1];
end
clr2 = smooth_colormap(clr2,5);
cmp = [clr1;clr2];
for ik=1:10
  cmp = smooth_colormap(cmp,25);
end

POS = [0.05 0.08 0.44 0.85; ...
       0.53 0.08  0.44 0.85];



LONN=[];
LATN=[];

%
% HYCOM:
% GoM region HYCOM: indices are for subsampled IAS HYCOM-TSIS domain!
GOM=[362   108
   288    50
   159     4
     9     5
     3   141
     8   401
   330   407
   398   405
   511   305
   520   210
   507   148
   429   118];

% GoM region, NEMO:
GOMN=[         100         365
         651         337
        1091         687
        1246         798
        1512         881
        1787         998
        1904        1292
        1710        1914
          23        1920
           8         748];




%Iexpt = find(EXON==1);
%for jj=1:length(Iexpt);
%  ixx = Iexpt(jj);
nmexp = EXPT(ixx).Name;
pthd1 = EXPT(ixx).path;
%  fmatout = sprintf('%shycom_LCcontour_%2.2i.mat',pthmat,ixx);
fmatout = sprintf('%shycom_t2Z%4.4i_hindcast%2.2i.mat',pthmat,abs(Z0),ixx);

fprintf('Loading %s\n',fmatout);
load(fmatout);
%
% Load saved LC contours corresponding to this T(z0)
fmatout2 = sprintf('%sHYCOM_dT%2.2i_Z%3.3i_contour_hind%2.2i.mat',...
                  pthmat,round(T0*10),abs(Z0),ixx);
fprintf('Loading %s\n',fmatout2);
load(fmatout2);

LCH=LCXY;
TM = TZH.TM;

% NEMO T contours:
fTNM = sprintf('%sNEMO_dT%2.2i_Z%3.3i_contour.mat',pthmat,round(T0*10),abs(Z0));
fprintf('Loading %s\n',fTNM);
load(fTNM);
LCN=LCXY;

figure('Position',[1300 529 1266 756]); clf;


dnmb = TM(irec);
dv = datevec(dnmb);
date_str = sprintf('%4.4i/%2.2i/%2.2i',dv(1:3));
fprintf('Plotting %s\n',datestr(dnmb));

FF = squeeze(TZH.Tz(irec,:,:));

icntr = find(LCH.TM==dnmb);
Xh=LCH.XY(icntr).X;
Yh=LCH.XY(icntr).Y;

if ~exist('INH','var');
  HHx = TZH.HH;
  [mm,nn]=size(HHx); 
  [XM,YM]=meshgrid([1:nn],[1:mm]);
  INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
  clear XM YM;

  LONH = TZH.LON;
  LATH = TZH.LAT;
end

% Subtract spatial mean T
dmm=FF;
dmm(INH==0)=nan;
tM=nanmean(nanmean(dmm));
thc = FF-tM;
thc(INH==0)=nan;



% 
% NEMO field:
% 
yr=dv(1);
mo=dv(2);
dm=dv(3);
dnmb1=datenum(yr,mo,1);
dnmb2=dnmb1+32;
v2=datevec(dnmb2);
dnmb2=datenum(v2(1),v2(2),1);
d2=dnmb2-datenum(yr,mo,1);

fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';
fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
fprintf('Reading NEMO: %s\n',fnemo);

fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);

if isempty(LONN)
  fmesh=sprintf('%smesh_mask.nc',fpwd);
  dmm = ncread(fmesh,'nav_lon');
  LONN = dmm';
  dmm = squeeze(ncread(fmesh,'nav_lat'));
  LATN = dmm';

  [mm,nn] = size(LONN);

  [XM,YM]=meshgrid([1:nn],[1:mm]);
  INN = inpolygon(XM,YM,GOMN(:,1),GOMN(:,2));
  clear XM YM

  ZZ = ncread(fin,'deptht');
  ZZ = -ZZ;

  dZ = abs(abs(ZZ)-abs(Z0));
  iz0 = find(dZ==min(dZ));
end

tnm = squeeze(ncread(fin,'toce',[1 1 iz0 dm],[Inf Inf 1 1]));
tnm = tnm';
I=find(tnm==0);
tnm(I)=nan;

% Subtract spatial mean ssh
dmm=tnm;
dmm(INN==0)=nan;
tM=nanmean(nanmean(dmm));
tnm = tnm-tM;
tnm(INN==0)=nan;





fnb=1;
c1=-3;
c2=3;
hhN=[];
xl1=-98;
xl2=-81;
yl1=18;
yl2=30.;
clb=1;

% ==============
% NEMO
% ==============
ps1=POS(1,:);
stt=sprintf('NEMO T %im, %s',abs(Z0),date_str);
sub_plot_2fld(ps1,tnm,LONN,LATN,hhN,c1,c2,stt,INN,xl1,xl2,yl1,yl2,clb);
contour(LONH,LATH,HHx,[0 0],'k');

Xn=LCN.XY(irec).X;
Yn=LCN.XY(irec).Y;
% Clean contours near Cuba:
xc1=-86.8;
xc2=-84.07;
xc3=xc2;
xc4=-81.1;
yc1=22.8;
yc2=yc1;
yc3=25.1;
yc4=yc3;
I=find(Yn<=yc1 & Xn>xc1);
Xn(I)=nan;
Yn(I)=nan;
I=find(Yn<=yc3 & Xn>xc3);
Xn(I)=nan;
Yn(I)=nan;
plot(Xn,Yn,'.');

% ==============
% HYCOM-TSIS
% ==============
ps1=POS(2,:);
stt2=sprintf('HYCOM-TSIS Hcst%i, T %im, %s',ixx,abs(Z0),date_str);
clb=0;

sub_plot_2fld(ps1,thc,LONH,LATH,HHx,c1,c2,stt2,INH,xl1,xl2,yl1,yl2,clb,'cmp',cmp);
exnm = EXPT(ixx).Name;
text(-92,18.5,exnm,'Fontsize',12);

Xh=LCH.XY(irec).X;
Yh=LCH.XY(irec).Y;
I=find(Yh<=yc1 & Xh>xc1);
Xh(I)=nan;
Yh(I)=nan;
I=find(Yh<=yc3 & Xh>xc3);
Xh(I)=nan;
Yh(I)=nan;
plot(Xh,Yh,'.');



btx = 'plot_T2z_hycom_nemo.m'; 
bottom_text(btx,'pwd',1);
