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

fmatout = sprintf('%sgomu2ias_indx.mat',pthmat);


%Read 0.03 IAS HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
HH(isnan(HH))=100;

% Topo, grid for 501 GOMu:
pthtopoR = '/nexsan/archive/GOMu0.04_501/topo/';
ftopoR = sprintf('%sdepth_GOMu0.04_03i.nc',pthtopoR);
HHR  = -1*squeeze(nc_varget(ftopoR,'depth'));
HHR(isnan(HHR))=100;
latR = nc_varget(ftopoR,'Latitude');
lonR = nc_varget(ftopoR,'Longitude');
[mm,nn]=size(HHR);
HHR(isnan(HHR))=100;

I = find(lonR>180.);
lonR(I) = lonR(I)-360.;

xlim1=min(lonR);
xlim2=max(lonR);
ylim1=min(latR);
ylim2=max(latR);

AA = HH;
AA(LON>xlim2-0.1)=nan;
AA(LAT<ylim1+0.1)=nan;
AA(LAT>ylim2-0.1)=nan;
AA(HH>=0)=nan;
IN = find(~isnan(AA));
lIN=length(IN);

[LONR,LATR] = meshgrid(lonR,latR);

INDX.IAS_HYCOM_Iocn = IN; 
INDX.IAS_HYCOM_idm  = nh;
INDX.IAS_HYCOM_jdm  = mh;


tic;
for ipp = 1:lIN
  if mod(ipp,10000)==0
    fprintf('  %4.1f%% done, %7.4f min ...\n',ipp/lIN*100,toc/60);
    tic;
  end

  ih=IN(ipp);
  x0=LON(ih);
  y0=LAT(ih);
  dx=abs(x0-lonR);
  dy=abs(y0-latR);
  i0=find(dx==min(dx),1);
  j0=find(dy==min(dy),1);
% Find 3 other points to close the grid bounding hycom point:
    dhx=lonR(i0)-x0;
  if abs(dhx)>0
      im1=i0-sign(dhx);
  else  % hycom lon=nemo lon
    im1=i0-1;
  end
  dhy=latR(j0)-y0;
  if abs(dhy)>0
      jm1=j0-sign(dhy);
  else
    jm1=j0-1;
  end

  vx=[lonR(i0),lonR(i0),lonR(im1),lonR(im1)];
  vy=[latR(j0),latR(jm1),latR(jm1),latR(j0)];
  inp=inpolygon(x0,y0,vx,vy);
  if ~inp,
    error('Check points, not inside NEMO grid ...\n');
  end

  f_chck=0;
  if f_chck==1
    clf;
    hold on;
    plot(LONR(j0-2:j0+2,i0-2:i0+2),LATR(j0-2:j0+2,i0-2:i0+2),'k.');


    plot(x0,y0,'g*');
    plot(LONR(j0,i0),LATR(j0,i0),'ro');
    plot(LONR(jm1,im1),LATR(jm1,im1),'bo');
    plot(vx,vy,'-');

   keyboard
  end

  INDX.I_GOMu(ipp,:)=[i0,i0,im1,im1];
  INDX.J_GOMu(ipp,:)=[j0,jm1,jm1,j0];

end

fprintf('Saving %s\n',fmatout);
save(fmatout,'INDX');





