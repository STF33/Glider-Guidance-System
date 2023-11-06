% Plot SSH NEMO 
%

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


%----------------------------------------
%
% Interpolate NEMO SSH onto HYCOM
%  
%----------------------------------------
fprintf('SSH interpolated into HYCOM-TSIS %s\n',datestr(dnmb));

% (1) interpolate NEMO to HYCOM
sshNi = sub_nemo2tsis_ssh(dnmb,pthnemo,pthtopo,pthdata,LAT,LON,HH,DX,DY);

% (2) interpolate GLORYS to HYCOM
%sshGi = sub_glorys2tsis_ssh(dnmb,pthnemo,pthtopo,pthdata,pthglorys,LAT,LON,HH);

% 
% Merge two fields:
%fnemo_indx = sprintf('%sNEMO_TO_HYCOM_TSIS_interp_pnts01.mat',pthdata);
%fprintf('Loading %s\n',fnemo_indx);

% Find NEMO bounding indices:
%NMI = load(fnemo_indx);
%Inm = NMI.IndxHYCOM;   % indices covered by NEMO
%[JH,IH] = ind2sub(size(HH),Inm);
%jbot = min(JH);
%jtop = max(JH);
%ilft = min(IH);
%irht = max(IH);

%Lmsk       = HH*0;
%Lmsk(HH<0) = 1;
%Lmsk(Inm)  = 0;

%sshN = sshGi;
%sshN(Inm) = sshNi(Inm);

sshN = sshNi;

% Smoothing along the NEMO OB is skipped - see interp_nemo/interp2D_nemo_glorys_hycom.m
%
% Demean:
I=find(INH==1);
sshM=nanmean(sshN(I));
sshN=sshN-sshM;

figure(1); clf;
set(gcf,'Position',[1573         560         949         751]);
axes('Position',[0.08 0.4 0.8 0.5]);
 % hindcast
pcolor(LON,LAT,sshN); shading flat;
colormap(cmp);
clm1=-0.3;
clm2=0.6;
%  caxis([-0.5 0.5]);
caxis([clm1, clm2]);
hold;
contour(LON,LAT,sshN,[0:0.2:1],'k-','linewidth',1);
contour(LON,LAT,sshN,[0.17 0.17],'r-','linewidth',1);


contour(LON,LAT,HH,[0 0],'-','Color',[0.5, 0.5, 0.5]);

axis('equal');
set(gca,'xlim',[-97.8 -80.7],...
 'ylim',[18.5 31],...
 'xtick',[-98:2:-82],...
 'ytick',[18:2:32],...
 'tickdir','out',...
 'Fontsize',12);


clb=colorbar('SouthOutside');
set(clb,'Position',[0.18 0.1 0.66 0.025],...
        'Fontsize',13,...
        'Ticks',[-1:0.1:1],...
        'Ticklength',0.025);
stl = sprintf('SSH NEMO %s',datestr(dnmb));
title(stl);

btx = 'plot_SSH_NEMO.m';
bottom_text(btx,'pwd',1);




















