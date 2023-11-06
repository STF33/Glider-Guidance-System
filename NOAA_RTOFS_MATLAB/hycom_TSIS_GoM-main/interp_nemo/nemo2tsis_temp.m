% INterpolate 1km NEMO onto HYCOM-TSIS
% Outside NEMO domain - use GLORYS
% Temperature

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
pthnemo = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
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
% First, interpolate NEMO to HYCOM
% Second GLORYS to HYCOM over domain outside NEMO
% Thrid, merge two interpolated fields
%
% Temp 
fldnm = 'toce';  % soce - S
nlr = 75; 
for iz=1:nlr
  TT = sub_get_NEMO_TS(dnmb,fldnm,iz);

  TT = sub_fill_land(TT);
  TTi = sub_interp_nemo2hycom(LONN,LATN,LON,LAT,HH,DX,DY,TT);

  TNEMO(iz).T=TTi;
		if f_mat>0
				foutp = sprintf('%sssh_nemo2hycom_%4.4i%2.2i%2.2i.mat',pthnemo,DV(1:3));
				fprintf('Saving %s\n',foutp);
				save(foutp,'TNEMO');
		end

end
  

% Get Glorys field:
[fglnm,flglrs] = sub_find_GLORYS_file(pthglorys,dnmb);

% Interpolate onto HYCOM-TSIS


