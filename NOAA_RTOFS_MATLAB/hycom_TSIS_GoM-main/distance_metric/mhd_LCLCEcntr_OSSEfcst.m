% Calculate MHD for the LC and LCE contours 
% OSSE forecast experiments
%
% Add/remove LCE based on best MHD
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
lonW = -89;  % LCEs west of lonW are discarded, < -98 - include all LCEs

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

%    Nlc1 = length(LCEN);
%    Nlc2 = length(LCE);  % max # of LCEs during the forecast at any given day
    Nlc1 = sub_numbLCE(LCEN,irc);
    Nlc2 = sub_numbLCE(LCE,ihc); % # of LCE on this date 
    [EComb, icomb] = sub_lce_comb(Nlc2);
    if Nlc2 ~= length(EComb)-1;
      error(' Check EComb - does not match Nlc2-1=%i',Nlc2-1);
    end

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

% For HYCOM forecast - try all possible combinations of LC/LCEs and find the best
% that matches the NEMO 
    for jlce=1:icomb   % adding eddies to the f/cast
      [Xhc,Yhc,flg] = sub_combineLCLCE_EComb(xh1,yh1,LCE,ihc,EComb,jlce);
      if jlce>1 & flg==0 % no LCE 
        MHD_LCLCE(jlce,1)=1.e20;
        continue;
      end

      P = [Xc,Yc];    % NEMO NR
      Q = [Xhc,Yhc];  % HYCOM fcst LC 
      mhd = modified_hausdorff_distance(P,Q,'geo');
      MHD_LCLCE(jlce,1) = mhd;        
    end
   
    mhd_fcst = min(min(MHD_LCLCE));
    [j0,i0] = find(MHD_LCLCE==mhd_fcst,1);
    MHD_LCLCE0 = MHD_LCLCE;
    mhd_fcst0 = mhd_fcst;
    MHD(ihc,1) = mhd_fcst;

%  Persistence contour
% Add/remove LCE for persistence:
    MHD_LCLCE = [];
    NlcP = sub_numbLCE(LCE,ihc0);
    [ECombP, icombP] = sub_lce_comb(NlcP);
    for jlce=1:icombP
      [Xprs,Yprs,flg] = sub_combineLCLCE_EComb(xPs,yPs,LCE,ihc0,ECombP,jlce);
      if jlce>1 & flg==0 % no LCE 
        MHD_LCLCE(jlce,1)=1.e20;
        continue;
      end

      P = [Xc,Yc];    % NEMO NR
      Q = [Xprs,Yprs];  % HYCOM persistence
      mhd = modified_hausdorff_distance(P,Q,'geo');
      MHD_LCLCE(jlce,1) = mhd;       
    end
    mhd_fcst = min(min(MHD_LCLCE));
    [j0,i0] = find(MHD_LCLCE==mhd_fcst,1);

    MHDprst(ihc,1) = mhd_fcst;


    itime0 = FRUNS(ixx).TimePeriod;
    irun0  = FRUNS(ixx).Nfcst;

    btx = 'mhd_LCLCEcntr_OSSEfcst.m';
    dnmb = TM(ihc);
    f_chck = 0;
    if f_chck == ihc & irun0 == 3 & itime0 == 1
      fnmb=10;
      sub_check_ssh_v2(xn,yn,xh1,yh1,LCEN,LCE,MHD_LCLCE0,icomb,ihc,irc,...
                    flg_lce1,HH,LON,LAT,fnmb);  
      dstr = datestr(dnmb);
      esim = EXPT(iFcst).Name;
      stl = sprintf('%s %s fcst day=%i, MHD=%6.3f',esim,dstr,ihc,mhd_fcst0);
      title(stl);
      bottom_text(btx,'pwd',1);

      fnmb=11;
      sub_check_ssh_v2(xn,yn,xPs,yPs,LCEN,LCE,MHD_LCLCE,icombP,ihc0,irc,...
                    flg_lce1,HH,LON,LAT,fnmb);
      dstr = datestr(dnmb);
      esim = EXPT(iFcst).Name;
      stl = sprintf('%s %s PERSIST day=%i, MHD=%6.3f',esim,dstr,ihc,mhd_fcst);
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
%  if Nhnd>=10
%    fmat1 = sprintf('%sMHD_LCLCE_nemo_persist_hycomfcst%2.2i-%2.2i%2.2i.mat',...
%                  pthmat2,FRUNS(ixx).Nhndcst,FRUNS(ixx).TimePeriod,FRUNS(ixx).Nfcst);
%  end

if f_save == 1 
  fprintf('Saving %s\n',fmat1);
  save(fmat1,'MHD','TM');
else
  fprintf('  NOT SAVING !!!! \n');
end
 
end;
%end

fprintf('All done\n');


