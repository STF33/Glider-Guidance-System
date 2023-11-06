% Plot SSH difference bteween NEMO and HYCOM Predict
% forecasts (experiments #10 - 16)
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

fldnm = 'ssh';

% Choose f/cast and time
iFcst = 10;
itime = 2;
irun  = 1; 
fday  = 91; % forecast day 1-100 (or 90)


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

%ts times:
FCST  = sub_fcstPrdct_info(iFcst);
RUN   = FCST.TIME0(itime).RUN(irun);
pthd1 = RUN.pthbin;
TM    = RUN.TM;
YDAY  = RUN.jday;
nrc   = length(TM);



dnmb = TM(fday);  % interpolation date
DV = datevec(dnmb);


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


%-------------------------------------------------------------------------------- 
%
%  Get HYCOM f/cast SSH:
%-------------------------------------------------------------------------------- 
yr  = DV(1);
mo  = DV(2);
dm  = DV(3);
iday= dnmb-datenum(yr,1,1)+1;

% Day 1 - initial field - take from the hindcast:
fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
dmean=1;
ssh = sub_getSSH_hycom(fina,finb,INH,dmean);

%sshn = sub_getSSH_nemo(DV(ii,:),LONN,LATN,INN,dmean);



%----------------------------------------
%
% Interpolate NEMO SSH onto HYCOM
%  
%----------------------------------------
fprintf('SSH interpolated into HYCOM-TSIS %s\n',datestr(dnmb));

% (1) interpolate NEMO to HYCOM
sshNi = sub_nemo2tsis_ssh(dnmb,pthnemo,pthtopo,pthdata,LAT,LON,HH,DX,DY);

% (2) interpolate GLORYS to HYCOM
sshGi = sub_glorys2tsis_ssh(dnmb,pthnemo,pthtopo,pthdata,pthglorys,LAT,LON,HH);

% 
% Merge two fields:
fnemo_indx = sprintf('%sNEMO_TO_HYCOM_TSIS_interp_pnts01.mat',pthdata);
fprintf('Loading %s\n',fnemo_indx);

% Find NEMO bounding indices:
NMI = load(fnemo_indx);
Inm = NMI.IndxHYCOM;   % indices covered by NEMO
[JH,IH] = ind2sub(size(HH),Inm);
jbot = min(JH);
jtop = max(JH);
ilft = min(IH);
irht = max(IH);

%Lmsk       = HH*0;
%Lmsk(HH<0) = 1;
%Lmsk(Inm)  = 0;

sshN = sshGi;
sshN(Inm) = sshNi(Inm);

% Smoothing along the NEMO OB is skipped - see interp_nemo/interp2D_nemo_glorys_hycom.m
%
% Demean:
I=find(INH==1);
sshM=nanmean(sshN(I));
sshN=sshN-sshM;

%
c1 = -0.6;
c2 = 0.6;
clr1 = colormap_hot(200);
clr2 = colormap_cold(200);
cmp = [flipud(clr2);clr1];


% Diff SSH:
dltSSH = ssh-sshN;

%
% Find mean dltSSH over the GoM:
mdSSH = nanmean(abs(dltSSH(I)));
rmse = sqrt(nansum(dltSSH(I).^2)./length(I));


Lmsk = HH*0;
Lmsk(HH<0)=1;
lcmp=[0 0 0; 1 1 1];


figure(1); clf;
set(gcf,'Position',[1682 600 834 648]);
axes('Position',[0.09 0.22 0.86 0.7]);
hold on;
%pcolor(LON,LAT,Lmsk); shading flat;
%colormap(lcmp);
%caxis([0 1]);
%freezeColors;

axis('equal');
set(gca,'tickdir','out',...
    'xlim',[-98 -80.1],...
    'xtick',[-100:2:-70],...
    'ylim',[18 30.7],...
    'ytick',[18:34],...
    'Color',[0.8 0.8 0.8],...
    'Fontsize',12);

pcolor(LON,LAT,dltSSH); shading flat;
caxis([c1 c2]);
colormap(cmp);

contour(LON,LAT,HH,[-200 -200],'k-','Color',[0.6 0.6 0.6],...
        'Linewidth',1.6);

clb=colorbar('SouthOutside');
set(clb,'Position',[0.18 0.1 0.66 0.025],...
        'Fontsize',13,...
        'Ticks',[-1:0.2:1],...
        'Ticklength',0.025);
stl = sprintf('SSH(HYCOM)-SSH(NEMO), F/cst %i, F/day=%i, %s',iFcst, fday,datestr(dnmb));
title(stl);

btx = 'plot_dltSSHPrdct_nemo_hycom.m';
bottom_text(btx,'pwd',1);




















