% Calculate MHD for the LC and LCE contours 
% in the Eastern GoM
% East of 90W
%
% for the HYCOM forecasts and NEMO
% for reference - persistence: day 0 of the free run
% is used for 2011 only
% extracted in hycom_TSIS/distance_metric/extr_lc_hycom_fcst.m
%
%  Forecasts: decision tree:
%
% H/cast #2 - Full 2D SSH                         Fcst 
%        #3 - AVISO SSH tracks only               Fcst 
%        #6 - 30th pnt T/S profiles GoM           not finished
%        #7 - AVISO + UGOS PIES (T/S profiles)    Fcst 
%        #8 - AVISO + extended PIES               Fcst 
%
%  Predictability experiments: interpolated NEMO+GLORYS --> HYCOM
%        #10 

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
iFcst=8; % #3, #7, #8 = match the hindcast numbers used for iniialization

lnW = -90; % cutoff longitude, allow max 1 LCE
LCE1 = 1; % when cobine controus, pick only 1 eastern most LCE

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';

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

for ixx=1:Mruns
  nmexp = FRUNS(ixx).Hind_name;
  Nhnd  = FRUNS(ixx).Nhndcst;

  flnm = FRUNS(ixx).matname;
  fmatout = sprintf('%s%s',pthmat,flnm);
  if Nhnd>=10
    fmatout = sprintf('%s%s',pthmat2,flnm);
  end
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

  fprintf('Calculating MHD %s\n',nmexp);
  for ihc=1:nrc
    if mod(ihc,30)==0,
      fprintf('  %6.2f%% done ...\n',ihc/nrc*100);
    end

    dm1= TM(ihc);
    irc= find(TMN==dm1);

    xn = LCN(1).XY(irc).X; % NEMO LC contour
    yn = LCN(1).XY(irc).Y;

    [Xc,Yc] = sub_combineLCLCE(xn,yn,LCEN,lnW,irc,LCE1);

% Save 1st contour NEMO - persistence reference  
    if ihc==1
      xPs=Xc;
      yPs=Yc;
    end
%
% Add/remove LCE based on MHD 
% if adding LCE imporves MHD - add, otherwise not
% same for deleting LCE west of lnW 
% the reason is that if LCE is very close to lnW this creates large errors
% if NEMO is slightly outside the lnW   
    if ~isempty(irc);
      xh1 = LCXY.XY(ihc).X;
      yh1 = LCXY.XY(ihc).Y;
%
% Force the LCE
     [Xhc,Yhc] = sub_combineLCLCE(xh1,yh1,LCE,-100,ihc,LCE1);

% Save HYCOM 1st contour for persistence
      if ihc==1
        xhPs = Xhc;
        yhPs = Yhc;
      end

      P = [Xc,Yc];    % NEMO
      Q = [Xhc,Yhc];  % HYCOM fcst LC + 1 LCE
      mhd1 = modified_hausdorff_distance(P,Q,'geo');

      Q0 = [xh1,yh1];  % HYCOM fcst - LC only
      mhd0 = modified_hausdorff_distance(P,Q0,'geo');

      if mhd0<mhd1,  % discard LCE, leave LC only
        mhd1=mhd0;
      end;

% Persistence using NEMO
      Q = [xPs,yPs];
      mhdPS = modified_hausdorff_distance(P,Q,'geo');

% Persistence - 1st day of HYCOM - not identically same as NEMO
      Q = [xhPs,yhPs];
      mhdPS_hycom = modified_hausdorff_distance(P,Q,'geo');

      if mhd1>100 & ihc<10
        fprintf('mhd1 %4.1f>100 ixx=%i, ihc=%i %s\n',mhd1,ixx,ihc,datestr(dm1));
        keyboard
      end

    else
      mhd1 = nan;
      mhdPS = nan;
      mhdPS_hycom = nan;
    end;

    MHD(ihc,1)  = mhd1;         % mhd fcst - nemo
    MHD(ihc,2)  = mhdPS;        % persistence nemo - nemo
    MHD(ihc,3)  = mhdPS_hycom;  % persistence hycom - nemo
    
%keyboard
  end
%keyboard    
  fprintf('  min/max MHD=%8.2f  %8.2f, min/max prst MHD=%8.2f %8.2f\n',...
          min(MHD(:,1)),max(MHD(:,1)),min(MHD(:,2)),max(MHD(:,2)));

  fmat1 = sprintf('%sMHD_LCLCE_nemo_persist_hycomfcst%2.2i-%2.2i%2.2i.mat',...
                  pthmat,FRUNS(ixx).Nhndcst,FRUNS(ixx).TimePeriod,FRUNS(ixx).Nfcst);
  if Nhnd>=10
    fmat1 = sprintf('%sMHD_LCLCE_nemo_persist_hycomfcst%2.2i-%2.2i%2.2i.mat',...
                  pthmat2,FRUNS(ixx).Nhndcst,FRUNS(ixx).TimePeriod,FRUNS(ixx).Nfcst);
  end
  
  fprintf('Saving %s\n',fmat1);
  save(fmat1,'MHD','TM');
%keyboard 
end;
%end

fprintf('All done\n');


