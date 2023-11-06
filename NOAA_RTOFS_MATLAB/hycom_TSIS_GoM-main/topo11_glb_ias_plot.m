% Plot bathymetry
%
% Extract topo from Global reanalysis and grid
% similar to IAS0.03 grid/domain
%
% GOFS3.1 GLBb0.08  reanalysis uses T11
%
% Add UGOS PIES / or extended PIES

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close all

% Add UGOS PIES (=1) or UGOS extended PIES (=2)
pthosse = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/OSSE/';
f_pies = 1; 

% Get T07 Global:
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthglb = '/nexsan/hycom/GLBb0.08_53X/topo/';
frga   = sprintf('%sregional.grid.a',pthglb);
frgb   = sprintf('%sregional.grid.b',pthglb);
fdptha = sprintf('%sdepth_GLBb0.08_11.a',pthglb);
fdpthb = sprintf('%sdepth_GLBb0.08_11.b',pthglb);

fidRGb = fopen(frgb,'r');  % read I,J from regional.grid.b
aa  = fgetl(fidRGb);
dmm = aa(2:8);
IDM = str2num(dmm);
aa = fgetl(fidRGb);
dmm = aa(2:8);
JDM = str2num(dmm);
IJDM = IDM*JDM;
fclose(fidRGb);
npad=4096-mod(IJDM,4096);

fprintf('IDM=%i, JDM=%i\n',IDM,JDM);

% read lon/lat from GLBb regional grid file
fidRGa = fopen(frga,'r');
[plon,count] = fread(fidRGa,IJDM,'float32','ieee-be');
fseek(fidRGa,4*(npad+IJDM),-1);
[plat,count] = fread(fidRGa,IJDM,'float32','ieee-be');

disp('Reading lat/lon for GLBb0.08 ...')
plon=(reshape(plon,IDM,JDM))';
plat=(reshape(plat,IDM,JDM))';

fclose(fidRGa);

[mg,ng]=size(plon);

I=find(plon>180);
plon(I)=plon(I)-360;


% TSIS IAS grid:
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;



% Find TSIS IAS grid in GLBb grid:
D=distance_spheric_coord(LAT(1,1),LON(1,1),plat,plon);
[j1,i1]=find(D==min(min(D)));
fprintf('Closest point in GLBb is %8.4f m away from TSIS(1,1)\n',min(min(D)));
D=distance_spheric_coord(LAT(1,end),LON(1,end),plat,plon);
[j2,i2]=find(D==min(min(D)));
fprintf('Closest point in GLBb is %8.4f m away from TSIS(1,end)\n',min(min(D)));
D=distance_spheric_coord(LAT(end,end),LON(end,end),plat,plon);
[j3,i3]=find(D==min(min(D)));
fprintf('Closest point in GLBb is %8.4f m away from TSIS(end,end)\n',min(min(D)));
D=distance_spheric_coord(LAT(end,1),LON(end,1),plat,plon);
[j4,i4]=find(D==min(min(D)));
fprintf('Closest point in GLBb is %8.4f m away from TSIS(end,1)\n',min(min(D)));


% read bathymetry from regional.depth.a
depth_fid=fopen(fdptha,'r');

%dmm=fread(depth_fid,6,'float32','ieee-be');
%aa=fread(depth_fid,IJDM,'float32','ieee-be');
%dmm=fread(depth_fid,npad,'float32','ieee-be');
%fseek(depth_fid,6*4*(npad+IJDM),-1) % <--- not needed
[h,count]=fread(depth_fid,IJDM,'float32','ieee-be');
y=find(h>1e10);
h(y)=nan;
h=reshape(h,IDM,JDM)';

h=-h;
h(isnan(h))=100;
HG=h; % GLBb topo

fclose(depth_fid);

di=200;
Hglb=HG(j1-di:j3+di,i1-di:i2+di);
LONglb=plon(j1-di:j3+di,i1-di:i2+di);
LATglb=plat(j1-di:j3+di,i1-di:i2+di);

x1=plon(j1,i1);
x2=plon(j3,i2);
y1=plat(j1,i2);
y2=plat(j3,i2);

Hglb(Hglb>0)=nan;
HH(HH>0)=nan;

%cmp=colormap(jet(360));
cmp = flipud(colormap_cold(360));
YR=2011;
if f_pies==1
  finp=sprintf('%sugos_mooring_%i.mat',pthosse,YR);
  fprintf('Loading %s\n',finp);
  A=load(finp);
  plat=A.ugos_lat;
  plon=A.ugos_lon;
%  sttl = sprintf('All AVISO, PIES, analysis date=%s',...
%     datestr(dnmb0));
elseif f_pies==2
  finp=sprintf('%sextd_mooring_2011.mat',pthosse);
  A=load(finp);
  plat=A.extd_lat;
  plon=A.extd_lon;
%  sttl = sprintf('All AVISO, extd PIES, analysis date=%s',...
%     datestr(dnmb0));
end


c1=-9000;
c2=0;
figure(1); clf;
set(gcf,'Position',[820 500 1063 895]);
axes('Position',[0.08 0.3 0.8 0.6]);
pcolor(LONglb,LATglb,Hglb); shading flat;
axis('equal');
colormap(cmp);
caxis([c1 c2]);
hold on;
plot(x1,y1,'m.');
plot(x2,y2,'m.');
plot(x1,y2,'m.');
plot(x2,y1,'m.');

set(gca,'tickdir','out',...
	'xlim',[-105 -40],...
	'ylim',[0 40],...
	'Color',[0 0 0]);

hb=colorbar;
set(hb,'Position',[0.9 0.3 0.015 0.6],...
       'TickLength',0.03,...
       'Fontsize',12);
title('HYCOM GLB0.08');

btx='topo11_glb_ias_plot.m';
bottom_text(btx,'pwd',1);


%plot HYCOM TSIS
figure(2); clf;
set(gcf,'Position',[850 450 1063 895]);
axes('Position',[0.08 0.3 0.8 0.6]);
pcolor(LON,LAT,HH); shading flat;
colormap(cmp);
hold on
if f_pies>0
  plot(plon,plat,'.','Color',[1,0.5,0.],'markersize',8);
end
caxis([c1 c2]);


axis('equal');
set(gca,'tickdir','out',...
	'xlim',[min(min(LON)) max(max(LON))],...
	'ylim',[min(min(LAT)) max(max(LAT))],...
	'Color',[0 0 0]);

hb=colorbar;
set(hb,'Position',[0.9 0.3 0.015 0.6],...
       'TickLength',0.03,...
       'Fontsize',12);
title('0.03 HYCOM TSIS');

bottom_text(btx,'pwd',1);


