% Plot SSH NEMO 
% on NEMO grid

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

fldnm = 'ssh';
dnmb = datenum(2011,07,09);  
DV = datevec(dnmb);


pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthout    = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat    = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

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
[mm,nn] = size(LONN);
%[DXN,DYN]=sub_dx_dy(LONN,LATN);
[XM,YM]=meshgrid([1:nn],[1:mm]);
% GoM region, NEMO:
GOMN=[         100         365
         651         337
        1091         687
        1246         798
        1512         881
        1787         998
        1904        1292
        1710        1914
          23        1920
           8         748];
INN = inpolygon(XM,YM,GOMN(:,1),GOMN(:,2));
clear XM YM



HH0 = LATN*0;
SCT = sub_GoM_xsct(LONN,LATN,HH0);

% SSH:
dmean = 1;  % do not demean SSH
sshN = sub_getSSH_nemo(DV,LONN,LATN,INN,dmean);



%cmp=flipud(colormap_cold(360));
% Colormap
ncc=100;
cl1=[0,0.3,0.6];
cl2=[1,1,1];
clr1=mix_2colors(cl1,cl2,ncc);

%cl1=cl2;
%cl2=[0.6,1,0.8];
cl1=cl2;
cl2=[0.6,0.9,0.7];
clr2=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[1,0.8,0.2];
clr3=mix_2colors(cl1,cl2,ncc);
cmp = [clr1;clr2;clr3];
cmp = smooth_colormap(cmp,5);



figure(1); clf;
set(gcf,'Position',[1573         560         949         751]);
axes('Position',[0.08 0.4 0.8 0.5]);
 % hindcast
pcolor(LONN,LATN,sshN); shading flat;
colormap(cmp);
clm1=-0.3;
clm2=0.6;
%  caxis([-0.5 0.5]);
caxis([clm1, clm2]);
hold;

contour(LONN,LATN,sshN,[0.17 0.17],'k-','linewidth',1.6);

% Plot xsection line:
XX=SCT.long;
YY=SCT.latd;
plot(XX,YY,'r.');

% Plot Distance hash-lines
DST  = cumsum(SCT.dX);
xplt = [0:200:max(DST)+10];
for idx=1:length(xplt)
	L0=xplt(idx);
	dlt= abs(DST-L0);
	jd = find(dlt==min(dlt));
	xl0= XX(jd);
	yl0= YY(jd);
	if idx==1
		plot(xl0,yl0,'k.','Markersize',20);
	else
		plot(xl0,yl0,'r.','Markersize',20);
	end

end

btx='plot_ssh_NEMO_NEMOgrid.m';
axis('equal');
set(gca,'xlim',[-97.8 -80.7],...
 'ylim',[18. 31],...
 'xtick',[-98:2:-82],...
 'ytick',[18:2:32],...
 'tickdir','out',...
 'color',[0 0 0],...
 'Fontsize',12);
hb=colorbar;
set(hb,'Position',[0.81 0.35 0.015 0.5],...
 'Ticks',[-1.0:0.2:1.0],...
 'Fontsize',12);

stl=sprintf('ssh, NEMO  %s',datestr(dnmb));
title(stl);

bottom_text(btx,'pwd',1,'Position',[0.08 0.3 0.4 0.04]);























