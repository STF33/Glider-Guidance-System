% HYCOM Predictability experiments: NEMO+GLORYS interpolated for
% initial fields
%
% Add/remove LCE based for NEMO NR and Predictability expts to get best MHD
% Calculate MHD for the LC and LCE contours 
% in the Eastern GoM
% East of 90W
%
% for the HYCOM forecasts and NEMO
% for reference - persistence: day 0 of the free run
% is used for 2011 only
% extracted in hycom_TSIS/distance_metric/extr_lc_hycom_fcst.m
%
%
%  Predictability experiments: interpolated NEMO+GLORYS --> HYCOM
%        #10, ..., #16

%  Forecast (based on Hindcast) group (which is used to initialize f/cst): 7 or 8 currently
%      |
%      V
%  TimePeriod during which the f/casts are initialized (May-June 2011, Jan-Feb 2012)
%     |
%     V
%  Forecasts (90 day f/csts shifted 7 days - 7 fcsts)
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
iFcst=16; % Predictability 1a - forecasts 10 - 16 


% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';

FCST = sub_fcstPrdct_info(iFcst);

Nhnd = iFcst;   % initial cond from the hindcast
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

for ixx=1:Mruns
  nmexp = FRUNS(ixx).Hind_name;
  Nhnd  = FRUNS(ixx).Nhndcst;

  flnm = FRUNS(ixx).matname;
  fmatout = sprintf('%s%s',pthmat2,flnm);

  fprintf('Loading %s\n',fmatout);
  load(fmatout);

  TM  = LCXY.TM;
  nrc = length(TM);
% Persistence contours: initial state NEMO
  xPs     = [];
  yPs     = [];
  xhPs    = [];
  yhPs    = [];
  MHD     = [];
  MHDprst = [];
  MHDprst = [];
  INEMO_best = [];
  INEMO_Pbest = [];
  ihc0    = -100;

  fprintf('Calculating MHD %s\n',nmexp);
  for ihc=1:nrc
    if mod(ihc,30)==0,
      fprintf('  %6.2f%% done ...\n',ihc/nrc*100);
      fprintf('  min/max MHD=%8.2f  %8.2f, min/max prst MHD=%8.2f %8.2f\n',...
          min(MHD(:,1)),max(MHD(:,1)),min(MHDprst),max(MHDprst));
    end

    dm1= TM(ihc);
    irc= find(TMN==dm1);

    xn = LCN.XY(irc).X; % NEMO LC contour
    yn = LCN.XY(irc).Y;


% Save 1st contour NEMO - persistence reference  
    if ihc0<=0
      ihc0=ihc;
    end
    xPs = LCXY.XY(ihc0).X;
    yPs = LCXY.XY(ihc0).Y;
    xh1 = LCXY.XY(ihc).X;
    yh1 = LCXY.XY(ihc).Y;

    Nlc1 = sub_numbLCE(LCEN,irc);
    [ECombN, icombN] = sub_lce_comb(Nlc1);  % Possible LC/LCE comb for NEMO

    Nlc2 = sub_numbLCE(LCE,ihc); % # of LCE on this date 
    [EComb, icomb] = sub_lce_comb(Nlc2);

    if Nlc2 ~= length(EComb)-1;
      error(' Check EComb - does not match Nlc2-1=%i',Nlc2-1);
    end

    NlcP = sub_numbLCE(LCE,ihc0);
    [ECombP, icombP] = sub_lce_comb(NlcP);

    lonW = -92;  % LCEs west of lonW are discarded, < -98 - include all LCEs

% NEMO NR - all LCEs 
% Add/remove 1 LCE based on MHD 
% if adding LCE imporves MHD - add, otherwise not
    MHD_LCLCE = zeros(icombN,1)*nan;
    MHD_Prs = MHD_LCLCE;
    flg_lce1 = [];
    for ilce=1:icombN   % LC/LCE combination for NEMO NR
      [Xc,Yc,flg] = sub_combineLCLCE_EComb(xn,yn,LCEN,irc,ECombN,ilce,'lonW',lonW);

      if flg < 0
        MHD_LCLCE(ilce) = 1.e30;  % no LCEs comb - west of lonW
        MHD_Prs(ilce) = 1.e30;
        continue
      end

% For HYCOM forecast - try all possible combinations of LC/LCEs and find the best
% that matches NEMO 
      [Xhc,Yhc,mhd_osse] = sub_best_contour_mhd(LCE,LCXY,Xc,Yc,ihc);
      MHD_LCLCE(ilce) = mhd_osse;

%  Persistence contour
% Add/remove LCE for persistence:
      [Xprs,Yprs,mhd_prs] = sub_best_contour_mhd(LCE,LCXY,Xc,Yc,ihc0);
      MHD_Prs(ilce) = mhd_prs;

    end

    mhd_min = min(MHD_LCLCE);
    iNbest = find(MHD_LCLCE==mhd_min,1);
    MHD(ihc,1) = mhd_min;
    INEMO_best(ihc,1) = iNbest;

    mhd_min = min(MHD_Prs);
    iPbest = find(MHD_Prs==mhd_min,1);
    MHDprst(ihc,1) = mhd_min;
    INEMO_Pbest(ihc,1) = iPbest;


  end

  MHD(:,2) = MHDprst;

  fprintf('  min/max MHD=%8.2f  %8.2f, min/max prst MHD=%8.2f %8.2f\n',...
          min(MHD(:,1)),max(MHD(:,1)),min(MHD(:,2)),max(MHD(:,2)));

  fmat1 = sprintf('%sMHD_LCLCE_nemo_persist_OSSEfcst%2.2i-%2.2i%2.2i.mat',...
                  pthmat,FRUNS(ixx).Nhndcst,FRUNS(ixx).TimePeriod,FRUNS(ixx).Nfcst);


  fprintf('Saving %s\n',fmat1);
  save(fmat1,'MHD','TM','INEMO_best','INEMO_Pbest');
%keyboard 
end;
%end

fprintf('All done\n');


