% Plot GLORYS sUV fields
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


fid = [];

zz0 = -0.5;

dnmb = datenum(2011,5,1);  % interpolation date
DV = datevec(dnmb);

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';

% Get GLORYS field:
[fglnm,flglrs] = sub_find_GLORYS_file(pthglorys,dnmb);

LONN = nc_varget(flglrs,'longitude');
LATN = nc_varget(flglrs,'latitude');
[LNN,LTN] = meshgrid(LONN,LATN);

ZZN  = nc_varget(flglrs,'depth');
ZZN = -ZZN;

%
% HYCOM grid:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
HH(isnan(HH))=100;
ln1 = min(min(LON));
ln2 = max(max(LON));
lt1 = min(min(LAT));
lt2 = max(max(LAT));
ilm1 = max(find(LONN<=ln1));
ilm2 = max(find(LONN<=ln2));
jlm1 = max(find(LATN<=lt1));
jlm2 = max(find(LATN<=lt2));


%

nlon = length(LONN);
nlat = length(LATN);
ndpth= length(ZZN);

iz0 = max(find(ZZN>=zz0));
if isempty(iz0), iz0=1; end;

fprintf('Reading GLORYS Depth layer %i\n',iz0);
F = squeeze(nc_varget(flglrs,'uo',[0 iz0-1 0 0],[1 1 nlat nlon]));
U = F(jlm1:jlm2,ilm1:ilm2);
F = squeeze(nc_varget(flglrs,'vo',[0 iz0-1 0 0],[1 1 nlat nlon]));
V = F(jlm1:jlm2,ilm1:ilm2);

S = sqrt(U.^2+V.^2);


nint=200;
c1=0;
c2=1.2;
%CMP = create_colormap3v2(nint,c1,c2);
%cmp1 = CMP.colormap;
CMP = create_colormap8(nint,c1,c2);
cmp = CMP.colormap;
%

[mm,nn] = size(S);
figure(1);
set(gcf,'Position',[1636 629 880 703]);
pcolor(S); shading flat;
caxis([c1 c2]);
colormap(cmp);

axis('equal');
set(gca,'xlim',[1 nn],'ylim',[1 mm]);

stl=sprintf('GLORYS U/V zz0=%7.1fm, %s',zz0,datestr(dnmb));
title(stl,'Interpreter','none');


clb = colorbar;
set(clb,'Position',[0.92 0.1 0.02 0.8]);

btx = 'plot_GLORYS_UV.m';
bottom_text(btx,'pwd',1);










