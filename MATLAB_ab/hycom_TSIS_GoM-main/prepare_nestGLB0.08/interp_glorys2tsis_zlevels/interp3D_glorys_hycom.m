% Wrapper script
% interpolation onto NEMO vertical layers
%
% Interpolation of 3D fields - T, S, U
% Interpolate GLORYS onto HYCOM-TSIS
% 
% Note all interpolation weights, points are precomputed 
% NEMO - HYCOM takes long to locate points
% 
% If relax HYCOM fields need to be created: 
% write interpolated fields into climatology HYCOM *a. *b files:
% see: write_hycom_clim.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


fid = [];

f_write = 0;
fldnm = 'temp';
%fldnm = 'saln';
%fldnm = 'uvel';
%fldnm = 'vvel';

Nlrs = 75;     % depth layers in NEMO

dnmb = datenum(2011,1,15);  % interpolation date
DV = datevec(dnmb);

fprintf('Interpolating: %s %s\n\n',fldnm,datestr(dnmb));

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';

% Read HYCOM topo:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
HH(isnan(HH))=100;


[DX,DY]=sub_dx_dy(LON,LAT);
%[XM,YM]=meshgrid([1:nn],[1:mm]);

fprintf('%s interpolating into HYCOM-TSIS %s\n',fldnm,datestr(dnmb));

% Load NEMO grid:
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

AmskN = [];  % T/S mask to fill bottom cells
AmskG = []; 
for iz0 = 1:Nlrs
%keyboard
  tic;
  zz0 = ZZN(iz0);  % NEMO depth
  fprintf('  %s iz=%i Z=%6.2f\n',fldnm,iz0,zz0);

% interpolate GLORYS to HYCOM on NEMO depth layer
  [Aglrs, AmskG] = sub_glorys2tsis_3d(dnmb,pthnemo,pthdata,pthglorys,LAT,LON,HH,...
                             fldnm,zz0,AmskG);

	Ahycom = Aglrs;
	% Fixed record length 
	[mm,nn]= size(HH);
	JDM = mm;
	IDM = nn;
	IJDM=IDM*JDM;
	npad=4096-mod(IJDM,4096);
	toto=ones(npad,1);

	if f_write==1
	% Fill land
		F = sub_fill_land(Ahycom);  
		F = F';
		F = reshape(F,IJDM,1);

	if isempty(fid)  
		fout = sprintf('%s%s_glorys2hycom_%2.2i%2.2i%4.4i.dat',...
																pthoutp,fldnm,DV(3),DV(2),DV(1));
		fid  = fopen(fout,'w','ieee-be');
	end

	fwrite(fid,F,'float32','ieee-be');
	fwrite(fid,toto,'float32','ieee-be');
	fprintf('Written files: %s\n',fout);
			
	end

  fprintf(' --------  1 record %8.4f min\n\n',toc/60);
end;   % depth levels

fclose(fid);

%
% Update links to *a, and create *b climatology files that are needed
% for relax files (interpolated onto HYCOM hybrid layers T&S)
% This replaces the code write_hycom_clim.m
if f_write==1 & (strncmp(fldnm,'temp',4) | strncmp(fldnm,'saln',4))
  sgm = 'sig2';  % pressure reference

  flnma = sprintf('%s_zlev_hycom_%s_%2.2i%2.2i%4.4i.a',fldnm,sgm,DV(3),DV(2),DV(1));
  fouta = sprintf('%s%s',pthoutp,flnma); 
  ib = max(strfind(fout,'/'));
  fnm_new = fout(ib+1:end);
  fprintf('Creating link %s ---> %s\n',fnm_new,flnma);
  D1 = pwd;

  s1=sprintf('cd %s',pthoutp(1:end-1));
  s2=sprintf('rm -f %s',flnma);
  s3=sprintf('ln -sf %s %s',fnm_new,flnma);
  s4=sprintf('cd %s',D1);

%  s1='who';
%  [stat,cout] = system(s1);
  eval(s1);
  system(s2);
  system(s3);
  eval(s4);

% Create *b file
  [JDM,IDM] = size(HH);
  IJDM = IDM*JDM;
  sub_create_climb(pthoutp,fldnm,flnma,dnmb,IDM,JDM,ZZN);

end;











