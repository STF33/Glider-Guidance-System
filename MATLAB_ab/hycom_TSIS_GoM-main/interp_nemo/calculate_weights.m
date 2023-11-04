% Find weights - surrounding NEMO nodes
% for interpolation onto HYCOM
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 1; % save mat; =2 - load saved and finish missing dates
s_par = 1; % parallel session
% 

dnmb = datenum(2011,5,1);  % interpolation date
DV = datevec(dnmb);


% 
pthnemo = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthdata = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthtopo = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';

% Read HYCOM topo:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);

[XM,YM]=meshgrid([1:n],[1:m]);

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
% Find NEMO points for all HYCOM domain:
if s_par>0
  delete(gcp('nocreate'))
  if exist('spool','var'),
    delete(spool);
  end

  spool = parpool('local');

end

fprintf('Finding NEMO points for HYCOM grid ...\n');
[ah1,ah2]=size(HH);
[an1,an2]=size(LONN);
ln1=min(min(LONN));
ln2=max(max(LONN));
lt1=min(min(LATN));
lt2=max(max(LATN));

I=find(HH<0 & ...
      (LON>=ln1 & LON<=ln2) & ...
      (LAT>=lt1 & LAT<=lt2));
nI = length(I);

np = 12;
INEMO = zeros(nI,np); 
%tic;
%nsum = zeros(nI,1);
parfor ii=1:nI
  atot=ii;
  if mod(atot,100)==0,
    fprintf(' %6.2f%% processed ...\n',atot/nI*100);
  end
%  fprintf('ii=%i\n',ii);

  i0=I(ii);
  [jh,ih] = ind2sub(size(HH),i0);

  x0=LON(i0);
  y0=LAT(i0);

  dx=0.5*DX(i0);
  dy=0.5*DY(i0);
  R=min([dx,dy]);

  D=distance_spheric_coord(LATN,LONN,y0,x0);

  J0=find(D<R);
  nj=length(J0);
 
  if nj<np
    J0(nj+1:np)=0; 
  elseif nj>np
    J0=J0(1:np);
  end
  INEMO(ii,:)=J0;
end



if exist('spool','var'),
  delete(spool);
end


if f_mat==1
  fout = sprintf('%sNEMO_TO_HYCOM_TSIS_interp_pnts.mat',pthdata)
  fprintf('Saving %s\n',fout);
  IndxHYCOM = I;
  save(fout,'INEMO','IndxHYCOM');
end


