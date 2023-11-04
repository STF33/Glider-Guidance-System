% Calculate MHD for the LC and LCE contours 
% OSSE forecast experiments
%
% Add/remove LCE based for NEMO NR and OSSE to get best MHD
%
% for the HYCOM forecasts and NEMO
% for reference - persistence: day 0 of the free run
% is used for 2011 only
% extracted in hycom_TSIS/distance_metric/extr_lc_hycom_OSSEfcst.m
%
%  Forecasts: decision tree:
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

f_save = 1;
% 1 forecast at a time
iFcst=2; % #2, #3, #6, #7, #8 = match the hindcast numbers used for iniialization

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

for ixx=1:Mruns
  nmexp = FRUNS(ixx).Hind_name;
  Nhnd  = FRUNS(ixx).Nhndcst;

  flnm = FRUNS(ixx).matname;
%  fmatout = sprintf('%s%s',pthmat,flnm);
%  if Nhnd>=10
    fmatout = sprintf('%s%s',pthmat2,flnm);
%  end
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
  INEMO_best = [];
  INEMO_Pbest = [];
  ihc0= -100;

  fprintf('Calculating MHD %s\n',nmexp);
  for ihc=1:nrc
    if mod(ihc,30)==0,
      fprintf('  %6.2f%% done ...\n',ihc/nrc*100);
      fprintf('  min/max MHD=%8.2f  %8.2f, min/max prst MHD=%8.2f %8.2f\n',...
          min(MHD(:,1)),max(MHD(:,1)),min(MHDprst),max(MHDprst));
    end

    dm1= TM(ihc);
    irc= find(TMN==dm1);
%
% Control run = Nature Run from NEMO
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

    Nlc1 = sub_numbLCE(LCEN,irc);  % # of LCEs 
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
% Try all LC/LCE combinations for NEMO and HYCOM OSSE 
% to get best MHD
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
    iPbest = find(MHD_Prs==mhd_min);
    MHDprst(ihc,1) = mhd_min;
    INEMO_Pbest(ihc,1) = iPbest;
    


    itime0 = FRUNS(ixx).TimePeriod;
    irun0  = FRUNS(ixx).Nfcst;

    btx = 'mhd_LCLCEcntr_best_OSSEfcst.m';
    dnmb = TM(ihc);
    f_chck = 0;
    if f_chck == ihc & irun0 == 1 & itime0 == 1
      fnmb=10;
     mhd_osse =  sub_check_ssh_v3(xn,yn,xh1,yh1,LCEN,LCE,LCXY,iNbest,ihc,irc,...
                    ECombN,HH,LON,LAT,fnmb);  
      dstr = datestr(dnmb);
      esim = EXPT(iFcst).Name;
      stl = sprintf('%s %s fcst day=%i, MHD=%6.3f',esim,dstr,ihc,mhd_osse);
      title(stl);
      bottom_text(btx,'pwd',1);

      fnmb=11;
      mhd_prst = sub_check_ssh_v3(xn,yn,xPs,yPs,LCEN,LCE,LCXY,iPbest,ihc0,irc,...
                       ECombN,HH,LON,LAT,fnmb);
      dstr = datestr(dnmb);
      esim = EXPT(iFcst).Name;
      stl = sprintf('%s %s PERSIST day=%i, MHD=%6.3f',esim,dstr,ihc,mhd_prst);
      title(stl);
      bottom_text(btx,'pwd',1);

      keyboard
    end
      

%keyboard
  end
%keyboard 
  MHD(:,2) = MHDprst;
     
  fprintf('  min/max MHD=%8.2f  %8.2f, min/max prst MHD=%8.2f %8.2f\n',...
          min(MHD(:,1)),max(MHD(:,1)),min(MHD(:,2)),max(MHD(:,2)));

  fmat1 = sprintf('%sMHD_LCLCE_nemo_persist_OSSEfcst%2.2i-%2.2i%2.2i.mat',...
                  pthmat,FRUNS(ixx).Nhndcst,FRUNS(ixx).TimePeriod,FRUNS(ixx).Nfcst);

  if f_save == 1 
    fprintf('Saving %s\n',fmat1);
    save(fmat1,'MHD','TM','INEMO_best','INEMO_Pbest');
  else
    fprintf('  NOT SAVING !!!! \n');
  end
   
end;
%end

fprintf('All done\n');


