% Update version: calculates MHD for different combinations of 
% LC/LCEs in the f/cast and control run and picks the 
% best f/cast based on lowest MHD
% This solves the problems of "jumps" when LCE is shed
%
% Calculate MHD for the LC and LCE contours 
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
%        #6 - 30th pnt T/S profiles GoM           Fcst
%        #7 - AVISO + UGOS PIES (T/S profiles)    Fcst 
%        #8 - AVISO + extended PIES               Fcst 
% 
% Do not use for experiments >=10 - these are Predictability experiments 1a/1b
%  Predictability experiments: interpolated NEMO+GLORYS --> HYCOM
%        #10 
%
%  Forecast (based on Hindcast) group (which is used to initialize f/cst): 2, 3, 6, 7, 8
%      |
%      V
%  TimePeriod during which the f/casts are initialized (May-June 2011, Jan-Feb 2012)
%     |
%     V
%  runs (90 day f/csts shifted 7 days - 7 fcsts)
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
iFcst = 2; % #3, #7, #8 = match the hindcast numbers used for iniialization
irun1 = 7; % runs started 7days apart
irun2 = 7;
itime1= 1;
itime2= 2;

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';

FCST = sub_fcst_info(iFcst);
Nhnd = FCST.Nhind;  % initial cond from the  hindcast #
hnd_name = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields

% All hindcast experiments
load('hycom_tsis_expts.mat');

% Array of forecast runs
icc=0;
for it=itime1:itime2
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
LCN   = LCXY;
LCEN  = LCE;  % LCE contours
TMN   = LCN(1).TM;  % Nemo
Nlcec = length(LCEN);

clear LCXY LCE

for ixx=1:Mruns
  nmexp = FRUNS(ixx).Hind_name;
  Nhnd  = FRUNS(ixx).Nhndcst;

  flnm = FRUNS(ixx).matname;
  fmatout = sprintf('%s%s',pthmat,flnm);
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
    if isempty(irc); 
      fprintf('NEMO contour is missing %s\n\n',datestr(dm1));
      continue;
    end

    xn = LCN(1).XY(irc).X; % NEMO LC contour
    yn = LCN(1).XY(irc).Y;

% Compare f/cast w/wout LCEs to the control run
% start with no eddies in the f/cast and control run
% and then add eddies to find the best agreement
% - lowest MHD
%
%      [Xc,Yc] = sub_combineLCLCE_v2(xn,yn,LCEN,irc,LCE1);

    Nlcec = length(LCEN);
    Nlcef = length(LCE);
%
% Add/remove LCE based on MHD 
% if adding LCE imporves MHD - add, otherwise not
% same for deleting LCE west of lnW 
% the reason is that if LCE is very close to lnW this creates large errors
% if NEMO is slightly outside the lnW   
    xh1 = LCXY.XY(ihc).X;
    yh1 = LCXY.XY(ihc).Y;

    MHD_LCLCE = [];

    for ilce=1:Nlcec+1 % adding eddies to the control run contour
      Neddy = ilce-1;
      [Xc,Yc,flgC] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc,Neddy);

      for jlce=1:Nlcef+1   % adding eddies to the f/cast
        if ilce>1 & flgC==0
          MHD_LCLCE(jlce,ilce)= 1e20; % no LCEs for this time recrod
          continue;
        end
        Neddy = jlce-1;
        [Xhc,Yhc,flg] = sub_combineLCLCE_EddyN(xh1,yh1,LCE,ihc,Neddy);

        P = [Xc,Yc];    % HYCOMcontrol run
        Q = [Xhc,Yhc];  % HYCOM fcst LC 
        mhd = modified_hausdorff_distance(P,Q,'geo');
        MHD_LCLCE(jlce,ilce) = mhd;

      end
    end
%
% To reduce "jumps" when the LCE shedding timing is off
    mhd_fcst = min(min(MHD_LCLCE));
    [j0,i0] = find(MHD_LCLCE==mhd_fcst,1);
    
    if ihc>1 & abs(mhd_fcst-MHD(ihc-1,1))>20
      fprintf('Big dlt MHD %6.4f\n',abs(mhd_fcst-MHD(ihc-1,1)));
%      keyboard;
    end;
%
% For persistence based on the best score
% Persistence using day=t0 from HYCOM control run
    if ihc==1
      Neddy = 0;  % no eddies for persistence - only LC contour
      [xPs,yPs,flg] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc,Neddy);
    end

    Q  = [xPs,yPs];
    mhdPs = [];
    dmmPs = [];
    for ilce=1:Nlcec+1
      Neddy = ilce-1;
      [Xc,Yc,flg] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc,Neddy);
      P = [Xc,Yc];    % HYCOMcontrol run
      dmmPs(ilce) = modified_hausdorff_distance(P,Q,'geo');
    end
    mhdPS = min(dmmPs);

    MHD(ihc,1)  = mhd_fcst;         % mhd fcst - hycom control run
    MHD(ihc,2)  = mhdPS;            % persistence hycom (t=t0)
    MHD(ihc,3)  = j0-1;             % how many LCEs 
    MHD(ihc,4)  = i0-1;


%keyboard
  end
%  keyboard

  fprintf('  min/max MHD=%8.2f  %8.2f, min/max prst MHD=%8.2f %8.2f\n',...
          min(MHD(:,1)),max(MHD(:,1)),min(MHD(:,2)),max(MHD(:,2)));

  fmat1 = sprintf('%sMHD_LCLCE_nemo_hycomfcst%2.2i-%2.2i%2.2i-V2.mat',...
                  pthmat,FRUNS(ixx).Nhndcst,FRUNS(ixx).TimePeriod,FRUNS(ixx).Nfcst);
  
  fprintf('Saving %s\n',fmat1);
  save(fmat1,'MHD','TM');
%keyboard 
end;
%end

fprintf('All done\n');


