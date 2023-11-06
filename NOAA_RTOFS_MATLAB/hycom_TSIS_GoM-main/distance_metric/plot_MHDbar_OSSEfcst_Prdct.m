% Plot bar diagram for  MHD LC/LCE  OSSE forecasts
% and Predictability 1a (NEMO interpolated as initial fields) 
%
% Before analysis:
% extract LC/LCE contour
% compute MHD
% 
% Calculated MHD for the LC contours in 
% !!!!  mhd_LCLCEcntr_nemo_fcsthycom.m <----- Old version
% New vesion: mhd_LCLCEcntr_OSSEfcst.m
%
% LC is extracted in hycom_TSIS/extr_lc_hycom_nemoV1.m
% NEMO LC is extracted in hycom_TSIS/extr_lc_ssh_nemo.m
%
% Hindcasts used for the forecast runs:
% H/cast #2 - Full 2D SSH                         Fcst 3
%        #3 - AVISO SSH tracks only               Fcst 4
%        #6 - 3D T/S 1/30 NEMO
%        #7 - AVISO + UGOS PIES (T/S profiles)    Fcst 1
%        #8 - AVISO + extended PIES               Fcst 3
% 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear


pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';

Ntime = 2;
% Forecasts < 10 - OSSE f/casts
%                 each f/cast group: 7 f/casts with 1 week shift
%           10-16 - predictability with NEMO+GLORYS interpolated into HYCOM
%                   1 run initialized at day 0 then 1 week shifted
%                   1 f/cast = 1 run, each f/cast is 1 week shifted
%                   so 10-16 f/casts are similar to 1 f/cast group 
%                   This is 1a experiments
% 
%     Do not use 1b experiments - when 1a run is used +/- 1,2 days
%     to initialize f/casts with perturbed ICs around day 0 of 1a run
%     day 0 - starts May 8, then shifts 1 week for each 1a
%
%    So, here I use only 1a Predict runs 10-16 as 1 f/cast group
%
IFCST=[2, 3, 6, 7, 8, 10]; % forecast groups
Nfgr=length(IFCST); % # of forecasts groups 

pthmat = pthmat1;

% All hindcast experiments
load('hycom_tsis_expts.mat');

%
% Combine MHD by forecast groups
% for OSSE f/casts
ifc=0;
MHD = sub_combine_MHD_OSSEfcst(IFCST,pthmat);

% Combine all Predictability 1a experiments into 1 f/cast group #10
MHD = sub_combine_MHD_prdct1a(MHD,pthmat1);

% Weekly pools
WK=[1:7:92]; % weekly pools
nwk=length(WK)-1;

% Pool all runs for the same f/cast group
% Average by months:
for itime=1:2
  for ifc=1:Nfgr
    POOL(ifc).Time(itime).pm1=[];
    POOL(ifc).Time(itime).pm2=[];
    POOL(ifc).Time(itime).pm3=[];
    PRST(ifc).Time(itime).pm1=[];
    PRST(ifc).Time(itime).pm2=[];
    PRST(ifc).Time(itime).pm3=[];

    for iwk=1:nwk
      POOL(ifc).Time(itime).WK(iwk).pw=[];
      PRST(ifc).Time(itime).WK(iwk).pw=[];
    end

  end
end


for itime=1:Ntime
  for ifc=1:Nfgr
% Pool all runs - MHD
% groups by months
    dmm  = MHD(ifc).Time(itime).MHD;
    POOL = sub_pool_month_week(POOL,dmm,ifc,itime,WK);
 
% Persistence:
    dmm  = MHD(ifc).Time(itime).MHD_prst;
    PRST = sub_pool_month_week(PRST,dmm,ifc,itime,WK);
 
  end
end
%

% Hindcast colors 
% match with hindcast for
% comparison
CLR = [0.4 0.8 0.2;...
       0   0.6 0.9; ...
       0.5   1   0.7; ...
       1   0.4 0.5; ... 
       0.  1  0.3; ...
       1   0   0.9; ...
       0.8 0   0.4; ...
       1   0.8 0; ...
       0.8 0.5 0; ...
       0.7 0.6 0.4; ...
       0.5 0.2 1];


btx='plot_MHDbar_OSSEfcst_Prdct.m';


% Plot: all F/csts in 1 group for 3 30-day time periods
% coordinates of the Bars relative to X=# of the 30-day period, 1,2,3
nfg=40;
anls_nm=sprintf('MHD(km) LC/LCE months');
%sub_plotMHDprst_barsrt_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,PRST,anls_nm);
% Plot persistence as separate bar next to expt:
%sub_plotMHD_prstbar_barsrt_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,PRST,anls_nm);
% Plot perst as a line IQR with median dot
sub_plotMHD_prstline_barsrt_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,PRST,anls_nm);

bottom_text(btx,'pwd',1);

% Plot weekly bars
f_week = 0;
if f_week==1
	nfg=41;
	anls_nm=sprintf('MHD(m) LC/LCE weeks');
	sub_plotMHDprst_barsrt_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,PRST,WK,anls_nm);
	bottom_text(btx,'pwd',1);
end

% Plot time series of MHD scores
f_chck=0;
if f_chck==1
	ifc=4;
	it=1;
  mo=1; % month 
  if mo==1
		aa=POOL(ifc).Time(it).pm1;
		bb=PRST(ifc).Time(it).pm1;	
  elseif mo==2
		aa=POOL(ifc).Time(it).pm2;
		bb=PRST(ifc).Time(it).pm2;	
  else
		aa=POOL(ifc).Time(it).pm3;
		bb=PRST(ifc).Time(it).pm3;	
  end

	figure(2); clf;
	plot(aa);
	hold on;
	plot(bb);
	md1=median(aa);
	md2=median(bb);

	fnm = MHD(ifc).Fcst_OSSE;
	stl=sprintf('%s md=%3.1f prst md=%3.1f',fnm,md1,md2);
	title(stl);
  legend({'Fcst','Persist'});
end








