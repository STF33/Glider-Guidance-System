% Find indices of NEMO grid
% corresponding to HYCOM grid:
% 4 closest points around HYCOM grid point
% the 1st point is the closest NEMO point
% the last is farthest from HYCOM grid point
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 1; % save mat; =2 - load saved and finish missing dates
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

fmatout = sprintf('%snemo2hycom_indx.mat',pthmat);

fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';

dnmb=datenum(2011,6,1);
DV = datevec(dnmb);
yr=DV(1);
mo=DV(2);
dnmb1=datenum(yr,mo,1);
dnmb2=dnmb1+32;
v2=datevec(dnmb2);
dnmb2=datenum(v2(1),v2(2),1);
d2=dnmb2-datenum(yr,mo,1);

fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
fprintf('Reading NEMO: %s\n',fnemo);
fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);

fmesh=sprintf('%smesh_mask.nc',fpwd);
dmm = ncread(fmesh,'nav_lon');
LONN = double(dmm');
dmm = squeeze(ncread(fmesh,'nav_lat'));
LATN = double(dmm');

[mm,nn] = size(LONN);

% NEMO grid: LON/LAT is fixed in i, j directions
lonn=LONN(1,:)';
latn=LATN(:,1);


%Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
fprintf('Reading hycom %s\n',ftopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
HH(isnan(HH))=100;

% Find HYCOM points that are inside NEMO domain
Xmn=min(min(LONN));
Ymn=min(min(LATN));
Xmx=max(max(LONN));
Ymx=max(max(LATN));
XV=[Xmn,Xmn,Xmx,Xmx];
YV=[Ymn,Ymx,Ymx,Ymn];


fprintf('Searching NEMO indiced for HYCOM grid \n');
INP=inpolygon(LON,LAT,XV,YV);
IN=find(INP==1 & HH<0); % ocean points only inside NEMO domain
lIN=length(IN);

INDX.info='HYCOM ocean indices inside NEMO';
INDX.info2='4 vertices of NEMO grid boupnding HYCOM pnt, 1st- closest';
INDX.HYCOM_size=[mh,nh];
INDX.HYCOM_ocean=IN;

tic;
for ipp=1:lIN
  if mod(ipp,10000)==0
    fprintf('  %4.1f%% done, %7.4f min ...\n',ipp/lIN*100,toc/60);
    tic;
  end

  ih=IN(ipp);
  x0=LON(ih);
  y0=LAT(ih);
  dx=abs(x0-lonn);
  dy=abs(y0-latn);
  i0=find(dx==min(dx),1);
  j0=find(dy==min(dy),1);
% Find 3 other points to close the grid bounding hycom point:
		dhx=lonn(i0)-x0;
  if abs(dhx)>0
  		im1=i0-sign(dhx);
  else  % hycom lon=nemo lon
    im1=i0-1;
  end
  dhy=latn(j0)-y0;
  if abs(dhy)>0
		  jm1=j0-sign(dhy);
  else
    jm1=j0-1;
  end

  vx=[lonn(i0),lonn(i0),lonn(im1),lonn(im1)];
  vy=[latn(j0),latn(jm1),latn(jm1),latn(j0)];
  inp=inpolygon(x0,y0,vx,vy);
  if ~inp,
    error('Check points, not inside NEMO grid ...\n');
  end

  f_chck=0;
  if f_chck==1
		clf;
		hold on;
		plot(LONN(j0-2:j0+2,i0-2:i0+2),LATN(j0-2:j0+2,i0-2:i0+2),'k.');

		plot(x0,y0,'g*');
		plot(LONN(j0,i0),LATN(j0,i0),'ro');
  plot(LONN(jm1,im1),LATN(jm1,im1),'bo');
  plot(vx,vy,'-');
  end

  INDX.I_NEMO(ipp,:)=[i0,i0,im1,im1];
  INDX.J_NEMO(ipp,:)=[j0,jm1,jm1,j0];

end

fprintf('Saving %s\n',fmatout);
save(fmatout,'INDX');








