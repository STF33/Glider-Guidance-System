addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

dnmb = datenum(2011,5,1);  % interpolation date
zz0 = -0.5;

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthintrp  = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthmat    = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/datamat/';


% Load NEMO grid:
% Get NEMO grid:
f_get_grid=0;
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
fprintf('Loading NEMO grid %s\n',fgrd);
load(fgrd);

iz0 = max(find(ZZN>=zz0));
if isempty(iz0), iz0=1; end;
zz0 = ZZN(iz0);

S = sub_get_NEMO_TS(dnmb,'soce',iz0);
[mm,nn] = size(S);

c1=33;
c2=38;
CMP = create_colormap7(200,c1,c2);
cmp = CMP.colormap;
cnt = CMP.intervals;
nint= length(cnt);

%
figure(1); clf;
set(gcf,'Position',[1636 629 880 703]);
pcolor(S); shading flat;
caxis([c1 c2]);
colormap(cmp);

axis('equal');

set(gca,'xlim',[1 nn],...
        'ylim',[1 mm],...
        'color',[0 0 0]);

stl=sprintf('NEMO S, zz0=%7.1fm, %s',zz0,datestr(dnmb));
title(stl);

clb = colorbar;
set(clb,'Position',[0.92 0.1 0.02 0.8]);

btx = 'plot_NEMO_S.m';
bottom_text(btx,'pwd',1);


