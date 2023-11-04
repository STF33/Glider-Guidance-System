% HYCOM Predictability experiments: NEMO+GLORYS interpolated for
% Experiments 1b: 
%  
% For each f/cast: 10, ..., 16
% there are: control run = iRun =1
% and  4 perturbed simulations, iRun=2,..., 5
%
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
iFcst=16; % Prdct1b f/casts 10, ..., 16 
itime1=1;
itime2=2;
irun1=2;  % irun=2,...,5
irun2=5;


lnW = -100; % cutoff longitude, allow max 1 LCE
LCE1 = 1; % when combine controus, pick only 1 eastern most LCE

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';


% All hindcast experiments
%load('hycom_tsis_expts.mat');

FCST = sub_fcstPrdct_info(iFcst);



% Array of forecast runs
%icc=0;
%for it=itime1:itime2
%  for irun=irun1:irun2  % note: iRun=1 - control run
%    icc=icc+1;
%    FRUNS(icc).iFcst = iFcst;
%    FRUNS(icc).TimePeriod = it;
%    FRUNS(icc).iRun = irun;
%    FRUNS(icc).matname = sprintf('hycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',iFcst,it,irun);
%    FRUNS(icc).Hind_name = FCST.Hindcast_Name;
%  end
%end
%
%
%Mruns = icc;

% YC point:
x0  = -85;
y0  =  22;

%if f_mhd==1
% Array with control run
%  fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat); OLD
%fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat1);

%fprintf('Loading %s\n',fmatout);
%load(fmatout);
%LCN=LCXY;
%LCEN=LCE;  % LCE contours
%TMN = LCN(1).TM;  % Nemo

%keyboard

clear LCXY LCE


iTime = 0;
for itime=itime1:itime2
% Control run - HYCOM without perturbation   
  iTime = itime; 
  irun0 = 1; % control run => irun=1
  flnm = sprintf('hycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',iFcst,itime,irun0);    
  fmatcntr = sprintf('%s%s',pthmat2,flnm);
  fprintf('Loading Control run %s\n',fmatcntr);
  load(fmatcntr);
  TMN = LCXY.TM;
  LCN = LCXY;
  LCEN = LCE;

  for irun=irun1:irun2
    RUN0  = FCST.TIME0(itime).RUN(1); % control run
    RUN   = FCST.TIME0(itime).RUN(irun);
    pthd1 = RUN.pthbin;
    TM    = RUN.TM;
    YDAY  = RUN.jday;
    nrc   = length(TM);
    day_shift = RUN.dayT0_true-RUN.prdct_time0;
    TM    = TM+day_shift;
    DV    = datevec(TM);

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
    fmat1 = sprintf('%sMHD_LCLCE_Prdct1b_%s.mat',pthmat2,fcstname);
    nmexp = fcstname;
  

    fprintf('\n\n %s %s\n',fcstname,fmat1);
    fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',iFcst,itime,irun);
    fprintf(' Input data: %s\n',pthd1);
    fprintf(' Output saved: %s\n',fmat1);

    flnm = sprintf('hycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',iFcst,iTime,irun);
    fmatfcst = sprintf('%s%s',pthmat2,flnm);
    fprintf('Loading f/cast run %s\n',fmatfcst);
    load(fmatfcst);



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

     xn = LCN(1).XY(irc).X; % control LC contour
     yn = LCN(1).XY(irc).Y;
%
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
    if ~isempty(irc);
      xh1 = LCXY.XY(ihc).X;
      yh1 = LCXY.XY(ihc).Y;

      MHD_LCLCE = [];
% 
% Start with no eddies only LC
% then add 1, only 1 eddy + LC is tested each time, more than 1 is not needed
% since jumps are due only to  LCE closest to LC, others
% that far away do not cause problems
      for ilce=1:Nlcec+1   % adding eddies to the control run contour
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
%
% For persistence based on the best score
% Persistence using day=t0 from HYCOM control run
      if ihc==1
%          Neddy = 0;  % no eddies for persistence - only LC contour
        xPs0 = xn; 
        yPs0 = yn;
        ircPs = irc;
%          [xPs,yPs,flg] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc,Neddy);
      end

%        Q  = [xPs,yPs];
      mhdPs = [];
      dmmPs = [];
      MHD_LCLCE = [];
      for ilce=1:Nlcec+1
        Neddy = ilce-1;
        [Xc,Yc,flgC] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc,Neddy);

        for jlce=1:Nlcef+1
          if ilce>1 & flgC==0
            MHD_LCLCE(jlce,ilce)= 1e20; % no LCEs for this time recrod
            continue;
          end
          Neddy = jlce-1;
          [xPs,yPs,flg] = sub_combineLCLCE_EddyN(xPs0,yPs0,LCEN,ircPs,Neddy);

          Q  = [xPs,yPs];
          P = [Xc,Yc];    % HYCOMcontrol run
          mhd = modified_hausdorff_distance(P,Q,'geo');
          MHD_LCLCE(jlce,ilce) = mhd;
        end
      end
      mhdPs = min(min(MHD_LCLCE));
      [jp0,ip0] = find(MHD_LCLCE==mhdPs,1);

      if mhd_fcst>60 
        fprintf('WARNING: high  mhd %4.1f %s\n',mhd_fcst,datestr(dm1));

%          Neddy = j0-1;
%          [Xhc,Yhc,flg] = sub_combineLCLCE_EddyN(xh1,yh1,LCE,ihc,Neddy);
%
%          Neddy = i0-1;
%          [Xc,Yc,flg] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc,Neddy);
%
%          figure(10); clf;
%          plot(Xhc,Yhc,'r.');
%          hold on;
%          plot(Xc,Yc,'.');
%          axis('equal');
%          keyboard
      end

    else
      mhd_fcst = nan;
      mhdPs = nan;
    end;

    MHD(ihc,1)  = mhd_fcst;         % mhd fcst - hycom control run
    MHD(ihc,2)  = mhdPs;            % persistence hycom (t=t0)
    MHD(ihc,3)  = j0-1;             % how many LCEs 
    MHD(ihc,4)  = i0-1;
       
  end

  fprintf('  min/max MHD=%8.2f  %8.2f, min/max prst MHD=%8.2f %8.2f\n',...
    min(MHD(:,1)),max(MHD(:,1)),min(MHD(:,2)),max(MHD(:,2)));
%keyboard
  fprintf('Saving %s\n',fmat1);
  save(fmat1,'MHD','TM');

  end
end;
%end

fprintf('All done\n');


