% Plot U/V form NEMO-GLORYS interpolated 
% into HYCOM-TSIS grid
% saved in archive file (nest)
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

flnm = sprintf('archv.%4.4i_%3.3i_00.a',iyr,iday);
fina = sprintf('%sarchv.%4.4i_%3.3i_00.a',pthoutp,iyr,iday);
finb = sprintf('%sarchv.%4.4i_%3.3i_00.b',pthoutp,iyr,iday);

%
[ZM,ZZ] = sub_zz_zm(fina,finb,HH);
zz0 = nanmean(nanmean(ZM(lplt,:,:)));

%
% Get barotrop:
[F,n,m,l] = read_hycom(fina,finb,'u_btrop');
F = squeeze(F);
F(F>huge)=nan;
ubtrop = F;

[F,n,m,l] = read_hycom(fina,finb,'v_btrop');
F = squeeze(F);
F(F>huge)=nan;
vbtrop = F;

[F,n,m,l] = read_hycom(fina,finb,'u-vel.','r_layer',lplt);
F = squeeze(F);
F(F>huge)=nan;
ubcl = F;

[F,n,m,l] = read_hycom(fina,finb,'v-vel.','r_layer',lplt);
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
set(gca,'xlim',[1 n],'ylim',[1  m]);

stl=sprintf('interp NEMO+GLORYS U/V to HYCOM, zz0=%7.1fm, %s',zz0,datestr(dnmb));
title(stl,'Interpreter','none');


clb = colorbar;
set(clb,'Position',[0.92 0.1 0.02 0.8]);

btx = 'plot_HYCOMrelax_UV.m';
bottom_text(btx,'pwd',1);





