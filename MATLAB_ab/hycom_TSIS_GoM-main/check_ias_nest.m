% Plot bathymetry
% for the IAS nest file
% prepared in /home/ddmitry/codes/anls_mtlb_utils/hycom_TSIS/prepare_nest_py/
%   bath_nest.py
%

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close


% Get T07 Global:
pthtopo = '/Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03/';
fltopo = 'regional.depth-nest';
fdptha = sprintf('%sregional.depth-nest.a',pthtopo);
fdpthb = sprintf('%sregional.depth-nest.b',pthtopo);
fgrd   = sprintf('%sregional.grid',pthtopo);

GRD = read_grid_bath(fgrd,fdptha);

HH = GRD.Topo;
HH(HH>1.e20)=nan;
HH=-1*HH;

LON = GRD.PLON;
LAT = GRD.PLAT;


cmp = flipud(colormap_cold(360));

c1=-9000;
c2=0;
figure(1); clf;
axes('Position',[0.08 0.3 0.8 0.6]);
pcolor(LON,LAT,HH); shading flat;
axis('equal');
colormap(cmp);
caxis([c1 c2]);
hold on;

set(gca,'tickdir','out',...
	'xlim',[-98 -56.1],...
	'ylim',[10 32],...
	'Color',[0 0 0]);

hb=colorbar;
set(hb,'Position',[0.88 0.3 0.02 0.6]);
title('HYCOM IAS nest');



