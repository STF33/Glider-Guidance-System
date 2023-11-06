% Plot NEMO vertical sections
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


fid = [];

%fldnm = 'temp';
fldnm = 'saln';
%fldnm = 'uvel'; % there is no u/v in relax files
%fldnm = 'vvel';

nlrs = 30;     % depth layers in HYCOM TSIS

dnmb = datenum(2011,5,1);  % interpolation date
DV = datevec(dnmb);

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthrelax  = '/Net/kronos/ddmitry/hycom/TSIS/relax/';

% 
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

%
% Find section similar to HYCOM-TSIS
% make sure that the points are inside NEMO
sct = 'Crb2';
fsect = sprintf('%shycom_tsis_xsct%s.mat',pthoutp, sct);
load(fsect);
slon = SCT.Lon;
slat = SCT.Lat;
nxx = length(slon);
xlmin = min(min(LONN));
xlmax = max(max(LONN));
ylmin = min(min(LATN));
ylmax = max(max(LATN));
isc=0;
clear IJ
for ik=1:nxx
  x0=slon(ik);
  y0=slat(ik);
  if x0<xlmin | x0>xlmax | ...
     y0<ylmin | y0>ylmax
    continue;
  end

  di=abs(LONN(100,:)-x0);
  ii0=find(di==min(di));
  dj=abs(LATN(:,ii0)-y0);
  jj0=find(dj==min(dj));
  isc=isc+1;
  IJ(isc,1)=ii0;
  IJ(isc,2)=jj0;
end



fprintf('%s interpolating into HYCOM-TSIS %s\n',fldnm,datestr(dnmb));
Nlrs = 75;     % depth layers in NEMO
for iz0 = 1:Nlrs
%keyboard
  tic;
  zz0 = ZZN(iz0);  % NEMO depth
  fprintf('  %s iz=%i Z=%6.2f\n',fldnm,iz0,zz0);

		switch(fldnm)
			case('temp')
				fldnemo = 'toce';
			case('saln')
				fldnemo = 'soce';
			case('uvel')
				fldnemo = 'uoce';
			case('vvel')
				fldnemo = 'voce';
		end

		INN =[];
		dmean = 0;  % do not demean SSH
		switch(fldnm)
			case({'temp','saln'})
				AA = sub_get_NEMO_TS(dnmb,fldnemo,iz0);
			case({'uvel','vvel'})
				AA = sub_get_NEMO_UV(dnmb,fldnemo,iz0);
		end
		AA0 = AA;

end
