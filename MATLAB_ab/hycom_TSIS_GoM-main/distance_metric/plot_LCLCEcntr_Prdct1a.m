% Plot LC/LCE contours used for MHD 
% for Predictability 1a (NEMO interpolated) forecast
% for specified date, experiment, run
% and NEMO NR
%  
% To match with OSSE f/casts:
% H/cast #2 - Full 2D SSH                         Fcst 
%        #3 - AVISO SSH tracks only               Fcst 
%        #6 - 30th pnt T/S profiles GoM           not finished
%        #7 - AVISO + UGOS PIES (T/S profiles)    Fcst 
%        #8 - AVISO + extended PIES               Fcst 
%
% Specify dates to plot and model f/cast day to find
% corresponding Predict1a run
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

% 1 forecast at a time
iday0 = 90;  % f/cast days 1, 30, 60, 90  
dnmbF = datenum(2012,4,13); % date of the f/cast
dnmb0 = dnmbF-iday0+1;  % Start date of the f/cast
itime0= 2;  % Time period 1 or 2
lonW  = -89; % check mhd_LCLCEcntr_OSSEfcst.m to have the same

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';


% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;


% YC point:
x0  = -85;
y0  =  22;

%if f_mhd==1
% Array with NEMO LC
%  fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat); OLD
fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat1);
fprintf('Loading %s\n',fmatout);
load(fmatout);
LCN=LCXY;
LCEN=LCE;  % LCE contours
TMN = LCN(1).TM;  % Nemo


clear LCXY LCE

% NEMO record #:
irc= find(TMN==dnmbF);
%
% Control run = Nature Run from NEMO
xn = LCN.XY(irc).X; % NEMO LC contour
yn = LCN.XY(irc).Y;
Nlc1 = sub_numbLCE(LCEN,irc); % # number of the LCEs in NEMO

% Save 1st contour NEMO - persistence reference  
ihc0=1;


% NEMO NR - all LCEs 
MHD_LCLCE = [];
flg_lce1 = [];
for ilce=1:Nlc1+1   % LC/LCE combination for NEMO NR
  Neddy = ilce-1;  % Neddy=0 - only LC
  [Xc,Yc,flgC] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc,Neddy,'lonW',lonW);
  if ilce>1
    flg_lce1(ilce-1)=flgC;
  end
  xn = Xc;
  yn = Yc;
end

% 
% Get predictability 1a forecast
% Find corresponding Predictability run with interpolated NEMO
ihc = iday0;
irun_prdct = 1;
iFcstI = sub_find_predict_fcast(dnmb0,itime0);
[LCEpr,LCXYpr] = sub_get_predictLCLCE(iFcstI,itime0,irun_prdct);
[Xprd,Yprd,mhd_prd] = sub_best_contour_mhd(LCEpr,LCXYpr,Xc,Yc,ihc);
%[Xprd,Yprd,Xc,Yc,mhd_prd] = sub_best_contour_mhd_v2(LCEpr,LCXYpr,xn,yn,ihc,LCEN,irc,lonW);

% Persistence Predict1a:
[Xprs,Yprs,mhd_prst] = sub_best_contour_mhd(LCEpr,LCXYpr,Xc,Yc,ihc0);
%[Xprs,Yprs,Xc,Yc,mhd_prst] = sub_best_contour_mhd_v2(LCEpr,LCXYpr,xn,yn,ihc0,LCEN,irc,lonW);


btx = 'plot_LCLCEcntr_Prdct1a.m';
dnmb = dnmbF;


%chyc = [0, 0.4, 0.8];
cprs = [0.5, 0.5, 0.5];
cnem = [0.9, 0.4, 0];
cprd = [0.0,0.4,0.8];

Lmsk = HH*0;
Lmsk(HH<0)=1;
cmp=[0,0,0; 1,1,1];


fprintf('Plotting contours ...\n');

figure(1); clf;
set(gcf,'Position',[1545, 665, 808, 667]);
axes('Position',[0.15, 0.25, 0.7,0.7]);
pcolor(LON,LAT,Lmsk); shading flat;
colormap(cmp);
freezeColors;
hold on;
contour(LON,LAT,HH,[-5000:500:-1],'Color',[0.9 0.9 0.9]);

%plot(Xhc,Yhc,'.','Color',chyc);
plot(Xc,Yc,'.','Color',cnem);
plot(Xprs,Yprs,'.','Color',cprs);
plot(Xprd,Yprd,'.','Color',cprd);

axis('equal');
set(gca,'xlim',[-96 -82],...
        'ylim',[20 30]);

dstr = datestr(dnmb);
esim = sprintf('Predict1a run=%i',irun_prdct); 
stl  = sprintf('%s %s fcast day=%i',esim,dstr,ihc);
title(stl);
%bottom_text(btx,'pwd',1);

axes('Position',[0.2 0.08 0.3 0.11]);
hold on

dxx=0.4;
dyy=0.3;

x1 = 1;
y1 = 1;
plot([x1 x1+dxx],[y1 y1],'linewidth',3,'Color',cnem);
text(x1+1.2*dxx,y1,'NEMO NR');

y1 = y1+dyy;
plot([x1 x1+dxx],[y1 y1],'linewidth',3,'Color',cprs);
text(x1+1.2*dxx,y1,sprintf('Persist (%5.1f)',mhd_prst));

y1 = y1+dyy;
plot([x1 x1+dxx],[y1 y1],'linewidth',3,'Color',cprd);
text(x1+1.2*dxx,y1,sprintf('NEMO-Interp (%5.1f)',mhd_prd));

ymx = y1;

axis('equal');
set(gca,'xlim',[1 2],...
        'ylim',[0.9 ymx+dyy],...
        'visible','off');

bottom_text(btx,'pwd',1);

 

