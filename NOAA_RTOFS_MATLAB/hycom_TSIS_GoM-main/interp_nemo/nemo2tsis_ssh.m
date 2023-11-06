% INterpolate 1km NEMO onto HYCOM-TSIS
% Outside NEMO domain - use GLORYS
%
% Weights (or points) for interpolation are derived and saved in 
% calculate_weights.m

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 0; % save mat; =2 - load saved and finish missing dates
% 

dnmb = datenum(2011,5,1);  % interpolation date
DV = datevec(dnmb);


% 
pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

% Read HYCOM topo:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
ah1=mm;
ah2=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);

[XM,YM]=meshgrid([1:nn],[1:mm]);

% Get NEMO grid:
f_get_grid=0;
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
if f_get_grid==1
  [LONN,LATN,ZZN] = sub_get_NEMO_grid(dnmb);
  LONN = double(LONN);
  LATN = double(LATN);
  ZZN = double(ZZN);
  fprintf('Saving grid %s\n',fgrd);
  save(fgrd,'LONN','LATN','ZZN');
else
  fprintf('Loading NEMO grid %s\n',fgrd);
  load(fgrd);
end 
[mm,nn] = size(LONN);
[DXN,DYN]=sub_dx_dy(LONN,LATN);

%
% NEMO interpolated weights/indices:
fnemo_indx = sprintf('%sNEMO_TO_HYCOM_TSIS_interp_pnts01.mat',pthdata);
fprintf('Loading %s\n',fnemo_indx);
NMI = load(fnemo_indx);

%keyboard

% ------------------------------
% Check that HYCOM points match
% saved indices:
% Note the order of points is random
% Derived by parallel job
%
[ah1,ah2]=size(HH);
[an1,an2]=size(LONN);
ln1=min(min(LONN));
ln2=max(max(LONN));
lt1=min(min(LATN));
lt2=max(max(LATN));

I1=find(HH<0 & ...
      (LON>=ln1 & LON<=ln2) & ...
      (LAT>=lt1 & LAT<=lt2));


Ihycom = NMI.IndxHYCOM;
if length(Ihycom)<length(I1)
  error('Saved HYCOM indices mismatch requested HYCOM points');
end

if length(Ihycom)==length(I1)
  if sum(Ihycom)~=sum(I1)
    fprintf('WARNING: Saved HYCOM indices mismatch requested HYCOM points\n');
  else
    fprintf(' HYCOM Index check ok\n');
  end
end


%
% First, interpolate NEMO to HYCOM
% Second GLORYS to HYCOM over domain outside NEMO
% Thrid, merge two interpolated fields
%
% SSH:
INN =[];
dmean = 0;  % do not demean SSH
ssh = sub_getSSH_nemo(DV,LONN,LATN,INN,dmean);

%
% 
ssh0 = ssh;
ssh = sub_fill_land(ssh);
sshi = sub_interp_nemo2hycom(LONN,LATN,LON,LAT,HH,DX,DY,NMI,ssh);

if f_mat>0
		foutp = sprintf('%sssh_nemo2hycom_%4.4i%2.2i%2.2i.mat',pthnemo,DV(1:3));
		fprintf('Saving %s\n',foutp);
		save(foutp,'sshi');
end

f_chck=1;
if f_chck==1;
  figure(1); clf;
  axes('Position',[0.07 0.2 0.4 0.7]);
  pcolor(ssh); shading flat;
  caxis([-0.4 0.4]);
  colorbar('SouthOutside');
  axis('equal');
  set(gca,'xtick',[],...
          'ytick',[]);
  title('NEMO');

  axes('Position',[0.53 0.2 0.4 0.7]);
  pcolor(sshi); shading flat;
  caxis([-0.4 0.4]);
  colorbar('SouthOutside');
  axis('equal');
  set(gca,'xlim',[1 700],...
          'ylim',[200 830]);
  set(gca,'xtick',[],...
          'ytick',[]);
%          'xlim',[1 ah2],...
%          'ylim',[1 ah1]);
  title('HYCOM-TSIS');
 
end 



