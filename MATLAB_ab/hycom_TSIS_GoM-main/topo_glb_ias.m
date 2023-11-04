% Extract topo from Global reanalysis and grid
% similar to IAS0.03 grid/domain
%
% GLBb0.08 19.0/19.1 reanalysis uses T07

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close


% Get T07 Global:
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthglb = '/nexsan/archive/GLBb0.08_190/topo/';
frga   = sprintf('%sregional.grid.a',pthglb);
frgb   = sprintf('%sregional.grid.b',pthglb);
fdptha = sprintf('%sdepth_GLBb0.08_07.a',pthglb);
fdpthb = sprintf('%sdepth_GLBb0.08_07.b',pthglb);

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

fmat=sprintf('%sGLBb0.08_T07_TSIS_IAS.mat',pthmat);
Hglb=HG(j1:j3,i1:i2);
LONglb=plon(j1:j3,i1:i2);
LATglb=plat(j1:j3,i1:i2);

save(fmat,'Hglb','LONglb','LATglb');

