% Plot U/V form NEMO-GLORYS interpolated 
% into HYCOM-TSIS grid
% from hycom restart file created from
% archv archive file (nest)
% archv --> restart - using HYCOM (Wallcraft) code
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear
close all

dnmb = datenum(2011,5,1);  % interpolation date
lplt = 1;  % layer to plot: 1,..., 30

% HYCOM:
rg=9806;  % convert pressure to depth, m
huge=1e20;

DV = datevec(dnmb);

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthintrp  = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthrstrt  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/restart/';
pthmat    = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/datamat/';


% Read HYCOM topo:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
HH(isnan(HH))=100;
[mm,nn]=size(HH);
[JDM,IDM] = size(HH);
IJDM = IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);


iyr = DV(1);
iday = dnmb-datenum(iyr,1,1)+1;

flnm = sprintf('restart_%4.4i_%3.3i_00.a',iyr,iday);
fina = sprintf('%srestart_%4.4i_%3.3i_00.a',pthrstrt,iyr,iday);
finb = sprintf('%srestart_%4.4i_%3.3i_00.b',pthrstrt,iyr,iday);

%
%[ZM,ZZ] = sub_zz_zm(fina,finb,HH);
%zz0 = nanmean(nanmean(ZM(lplt,:,:)));

%
% Get barotrop:
[F,n,m,l] = read_hycom_restart(fina,finb,'ubavg',IDM,JDM);
F = squeeze(F);
F(F>huge)=nan;
ubtrop = F;

[F,n,m,l] = read_hycom_restart(fina,finb,'vbavg',IDM,JDM);
F = squeeze(F);
F(F>huge)=nan;
vbtrop = F;

[F,n,m,l] = read_hycom_restart(fina,finb,'u',IDM,JDM,'r_layer',lplt);
F = squeeze(F);
F(F>huge)=nan;
ubcl = F;

[F,n,m,l] = read_hycom_restart(fina,finb,'v',IDM,JDM,'r_layer',lplt);
F = squeeze(F);
F(F>huge)=nan;
vbcl = F;

U = ubcl+ubtrop;
V = vbcl+vbtrop;
S = sqrt(U.^2+V.^2);

c1=0.;
c2=1.2;
CMP = create_colormap8(200,c1,c2);
cmp = CMP.colormap;
cnt = CMP.intervals;
nint= length(cnt);

figure(1); clf;
set(gcf,'Position',[1636 629 880 703]);
pcolor(S); shading flat;
caxis([c1 c2]);
colormap(cmp);

axis('equal');
set(gca,'xlim',[1 nn],'ylim',[1  mm]);

stl=sprintf('restart NEMO+GLORYS U/V to HYCOM, ilr %i, %s, %s',...
             lplt,datestr(dnmb),flnm);
title(stl,'Interpreter','none');

clb = colorbar;
set(clb,'Position',[0.92 0.1 0.02 0.8]);

btx = 'plot_HYCOMrestart_UV.m';
bottom_text(btx,'pwd',1);




