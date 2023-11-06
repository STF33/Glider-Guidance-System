% Get demeaned ssh NEMO for the Gulf of Mexico
% SSH interpolated onto HYCOM IAS grid
% 
function sshN = sub_get_sshNEMO(dnmb);

fldnm = 'ssh';
%dnmb = datenum(2012,03,20);
DV = datevec(dnmb);

fprintf('Deriving demeaned SSH NEMO %s\n',datestr(dnmb));
fprintf(' SSH NEMO interpolated onto IAS HYCOM grid \n');

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthout    = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat    = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

% Info about All hindcast experiments
load('hycom_tsis_expts.mat');

% Load NEMO 2 hycom indices:
fmatindx = sprintf('%snemo2hycom_indx.mat',pthmat);
fprintf('Loading %s\n',fmatindx);
load(fmatindx);

% Read HYCOM topo:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
HH(isnan(HH))=100;
[mh,nh]=size(HH);
m=mh;
n=nh;

[DX,DY]=sub_dx_dy(LON,LAT);
%[XM,YM]=meshgrid([1:nn],[1:mm]);

% GoM region HYCOM:
GOM=[366   489
   476   531
   583   560
   576   646
   508   827
   336   848
   204   829
    64   798
    19   746
    16   662
    12   578
    25   455
    71   382
   165   356
   281   400];

[XM,YM]=meshgrid([1:n],[1:m]);
INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
Z0 = -10;
Iocn = find(INH==1 & HH<Z0);
Nocn=length(Iocn);
clear XM YM

%----------------------------------------
%
% Interpolate NEMO SSH onto HYCOM
%  
%----------------------------------------
fprintf('SSH interpolated into HYCOM-TSIS %s\n',datestr(dnmb));

% (1) interpolate NEMO to HYCOM
sshNi = sub_nemo2tsis_ssh(dnmb,pthnemo,pthtopo,pthdata,LAT,LON,HH,DX,DY);

sshN = sshNi;

% Smoothing along the NEMO OB is skipped - see interp_nemo/interp2D_nemo_glorys_hycom.m
%
% Demean:
I=find(INH==1);
sshM=nanmean(sshN(I));
sshN=sshN-sshM;

return





