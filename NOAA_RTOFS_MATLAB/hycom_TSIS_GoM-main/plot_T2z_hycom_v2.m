% Plot HYCOM interpolated T onto z level
% interpolation performed in: interp_hycom2z.m
% and extracted LC / LCE contours: extr_lc_temp_hycom.m
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
irec = 189;  % day

%EXON = zeros(9,1);
%EXON(9) = 1; % select expt to plot
% Select hindcast experiment #: 2 -- 9
ixx = 9;  % OSSE hindcast/free run # expt



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

POS = [0.1, 0.1, 0.8, 0.8];


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


figure('Position',[1300 529 800 800]); clf;


dnmb = TM(irec);
dv = datevec(dnmb);
date_str = sprintf('%4.4i/%2.2i/%2.2i',dv(1:3));
fprintf('Plotting %s\n',datestr(dnmb));

FF = squeeze(TZH.Tz(irec,:,:));

icntr = find(LCH.TM==dnmb);
Xh=LCH.XY(icntr).X;
Yh=LCH.XY(icntr).Y;

if ~exist('INH','var');
  HH = TZH.HH;
  [mm,nn]=size(HH); 
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



fnb=1;
c1=-3.5;
c2=3.5;
hhN=[];
xl1=-97.5;
xl2=-81.5;
yl1=18.6;
yl2=29.8;
clb=1;


Lmsk = HH*0;
Lmsk(HH<0)=1;

% ==============
% HYCOM-TSIS
% ==============
ps1=POS;
stt2=sprintf('HYCOM-TSIS Hcst%i, T %im, %s',ixx,abs(Z0),date_str);
clb=3;

sub_plot_dTz(ps1,thc,LONH,LATH,HH,c1,c2,stt2,INH,xl1,xl2,yl1,yl2,clb,Lmsk,'cmp',cmp);
exnm = EXPT(ixx).Name;
text(-92,32,exnm,'Fontsize',12);


Xh=LCH.XY(irec).X;
Yh=LCH.XY(irec).Y;
% Clean contours near Cuba:
f_cl=0;
if f_cl==1
	xc1=-86.8;
	xc2=-84.07;
	xc3=xc2;
	xc4=-81.1;
	yc1=22.8;
	yc2=yc1;
	yc3=25.1;
	yc4=yc3;

	I=find(Yh<=yc1 & Xh>xc1);
	Xh(I)=nan;
	Yh(I)=nan;
	I=find(Yh<=yc3 & Xh>xc3);
	Xh(I)=nan;
	Yh(I)=nan;
end
plot(Xh,Yh,'.');

btx = 'plot_T2z_hycom_v2.m'; 
bottom_text(btx,'pwd',1);
