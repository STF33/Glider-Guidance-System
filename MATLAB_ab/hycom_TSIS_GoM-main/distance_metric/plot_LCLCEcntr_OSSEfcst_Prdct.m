% Plot LC/LCE contours used for MHD 
% tunred off ---> with persistence and Predictability (NEMO interpolated) forecast
% see separate code to plot Predict1a plot_LCLCEcntr_Prdct1a.m
%
% for specified date, experiment, run
%
% H/cast #2 - Full 2D SSH                         Fcst 
%        #3 - AVISO SSH tracks only               Fcst 
%        #6 - 30th pnt T/S profiles GoM           not finished
%        #7 - AVISO + UGOS PIES (T/S profiles)    Fcst 
%        #8 - AVISO + extended PIES               Fcst 
%
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

% 1 forecast at a time
iFcst = 6; % #2, #3, #6, #7, #8 = match the hindcast numbers used for iniialization
lonW  = -86; % check mhd_LCLCEcntr_OSSEfcst.m to have the same
irun0 = 3; 
itime0= 1; % Time 1 = 2011, 2=2012
%iday0 = 90;  % days 1, 30, 60, 90  
iday0 = 90;

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


if iFcst<10
  FCST = sub_fcst_info(iFcst);
else
  FCST = sub_fcstPrdct_info(iFcst);
end
Nhnd = FCST.Nhind;  % initial cond from the  hindcast #
hnd_name = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
irun1 = FCST.run1;
irun2 = FCST.run2;
ntime=FCST.ntime; % 2 time windows for forecasts
NTime=ntime;

% All hindcast experiments
load('hycom_tsis_expts.mat');

% Array of forecast runs
icc=0;
for it=1:NTime
  for irun=irun1:irun2
    icc=icc+1;
    FRUNS(icc).Nhndcst = Nhnd;
    FRUNS(icc).TimePeriod = it;
    FRUNS(icc).Nfcst = irun;
    FRUNS(icc).matname = sprintf('hycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',Nhnd,it,irun);
    FRUNS(icc).Hind_name = EXPT(Nhnd).Name;
  end
end

%fmat = 'hycom_info_forecasts.mat';
%fprintf('saving forecast info %s\n',fmat);
%save(fmat,'FRUNS');

Mruns = icc;

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

% find time period and forecast run
ixx = [];
for ix0=1:Mruns
  if FRUNS(ix0).TimePeriod == itime0 & ...
     FRUNS(ix0).Nfcst == irun0 
    ixx = ix0;
    break;
  end
end



nmexp = FRUNS(ixx).Hind_name;
Nhnd  = FRUNS(ixx).Nhndcst;

flnm = FRUNS(ixx).matname;
fmatout = sprintf('%s%s',pthmat2,flnm);
fprintf('Loading %s\n',fmatout);
load(fmatout);

TM  = LCXY.TM;
nrc = length(TM);
% Persistence contours: initial state NEMO
xPs = [];
yPs = [];
xhPs= [];
yhPs= [];
MHD = [];
MHDprst = [];
ihc0= -100;

fprintf('Calculating MHD %s\n',nmexp);
ihc = iday0;

% 1st day of the f/cast:
dnmb0 = TM(1);

dm1= TM(ihc);
irc= find(TMN==dm1);
%
% Control run = Nature Run from NEMO
xn = LCN.XY(irc).X; % NEMO LC contour
yn = LCN.XY(irc).Y;

% Save 1st contour NEMO - persistence reference  
ihc0=1;

Nlc1 = sub_numbLCE(LCEN,irc);
[ECombN, icombN] = sub_lce_comb(Nlc1);  % Possible LC/LCE comb for NEMO

% NEMO NR - all LCEs 
% Add/remove 1 LCE based on MHD 
% if adding LCE imporves MHD - add, otherwise not
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


% Get best LC/LCE contour that results in lowest MHD for OSSE f/casts
[Xhc,Yhc,mhd_osse] = sub_best_contour_mhd(LCE,LCXY,Xc,Yc,ihc);
%[Xhc,Yhc,Xc,Yc,mhd_osse] = sub_best_contour_mhd_v2(LCE,LCXY,xn,yn,ihc,LCEN,irc,lonW);

% Get best LC/LCE contour that results in lowest MHD for persistence
[Xprs,Yprs,mhd_prst] = sub_best_contour_mhd(LCE,LCXY,Xc,Yc,ihc0);
%[Xprs,Yprs,Xc,Yc,mhd_prst] = sub_best_contour_mhd_v2(LCE,LCXY,xn,yn,ihc0,LCEN,irc,lonW);

% 
% Get predictability 1a forecast
% Find corresponding Predictability run with interpolated NEMO
%irun_prdct = 1;
%iFcstI = sub_find_predict_fcast(dnmb0,itime0);
%[LCEpr,LCXYpr] = sub_get_predictLCLCE(iFcstI,itime0,irun_prdct);
%[Xprd,Yprd,mhd_prd] = sub_best_contour_mhd(LCEpr,LCXYpr,Xc,Yc,ihc);


btx = 'plot_LCLCEcntr_OSSEfcst_Prdct.m';
dnmb = TM(ihc);


chyc = [0, 0.4, 0.8];
cprs = [0.5, 0.5, 0.5];
cnem = [0.9, 0.4, 0];
cprd = [0.0,0.9,0.4];

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

plot(Xhc,Yhc,'.','Color',chyc);
plot(Xc,Yc,'.','Color',cnem);
plot(Xprs,Yprs,'.','Color',cprs);  % persistence


axis('equal');
set(gca,'xlim',[-96 -82],...
        'ylim',[20 30]);


dstr = datestr(dnmb);
esim = EXPT(iFcst).Name;
stl = sprintf('%s %s fcast day=%i',esim,dstr,ihc);
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
plot([x1 x1+dxx],[y1 y1],'linewidth',3,'Color',chyc);
text(x1+1.2*dxx,y1,sprintf('OSSE fcast (%5.1f)',mhd_osse));

%y1 = y1+dyy;
%plot([x1 x1+dxx],[y1 y1],'linewidth',3,'Color',cprd);
%text(x1+1.2*dxx,y1,sprintf('NEMO-Interp (%5.1f)',mhd_prd));

ymx = y1;

axis('equal');
set(gca,'xlim',[1 1.8],...
        'ylim',[0.9 ymx+dyy],...
        'visible','off');

bottom_text(btx,'pwd',1);

 

