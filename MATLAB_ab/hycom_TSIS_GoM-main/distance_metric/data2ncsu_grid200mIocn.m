% Prepare data for Ruoying/ NCSU
% Expt 1a - NEMO-HYCOM 
% TOPO, LON/LAT, Iocn for 200m SSH RMSE
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


%Z0 = -10; % discard close to coastline points
Z0 = -200;

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/data2ncsu/'; 
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';


%Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
HH(isnan(HH))=100;

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
Iocn = find(INH==1 & HH<Z0);
Nocn=length(Iocn);
clear XM YM


GRID.Name        = 'HYCOM grid and ocean indices in the GoM inside z0 ';
GRID.cutoff_z0   = Z0;
GRID.topo        = HH;
GRID.LON         = LON;
GRID.LAT         = LAT;
GRID.Iocn        = Iocn;

flgrid = sprintf('%sHYCOM_grid_Iocn%4.4im.mat',pthmat,abs(Z0));
fprintf('Saving %s\n',flgrid);
save(flgrid,'GRID');


