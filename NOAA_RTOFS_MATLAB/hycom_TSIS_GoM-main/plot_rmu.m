% Plot relax fields for HYCOM-TSIS 
% nested within GLORYS
% to macth NEMO
% NEMO domain is smaller and relaxation out the NEMO domain is tronger
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close


pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthin = '/nexsan/people/abozec/TSIS/IASx0.03/topo/';
btx = 'plot_rmu.m';


% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

IDM = nn;
JDM = mm;
IJDM = nn*mm;

fina = sprintf('%srmu_nemo_strong.a',pthin);
finb = sprintf('%srmu_nemo_strong.b',pthin);

%fida = fopen(fina,'r','ieee-be');
%fidb = fopen(finb,'r');

% Read relaxaion file:
fida = fopen(fina,'r');
A = fread(fida,IJDM,'float32','ieee-be');
A = (reshape(A,IDM,JDM))';
fclose(fida);

A(A==0) = nan;
RL = 1./A*(1/3600); % sec^(-1) -> hours


figure(1); clf;
contour(LON,LAT,HH,[0 0],'k');
hold on;
pcolor(LON,LAT,RL); shading flat;
%  caxis([0  c2]);
axis('equal');

caxis([0 5]);

colorbar;
%Is=max(strfind(ftrca,'/'));
spp=sprintf('HYCOM-TSIS, Min Relax Time %6.2f hrs',min(min(RL)));
title(spp,'interpreter','none','Fontsize',12);

bottom_text(btx,'pwd',1,'Position',[0.08 0.1 0.4 0.05]);




