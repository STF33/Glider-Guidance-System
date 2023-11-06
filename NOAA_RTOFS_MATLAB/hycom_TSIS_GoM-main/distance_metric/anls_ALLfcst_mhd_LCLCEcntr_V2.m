% Analyze MHD for the forecasts
%
% Before analysis:
% extract LC/LCE contour
% compute MHD
% 
% Calculated MHD for the LC contours in 
% mhd_LCLCEcntr_nemo_fcsthycom.m
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
% Plot time series:
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

% Figure to plot:
f_mhd_tser = 0; % plot MHD time series 


pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
btx = 'anls_ALLfcst_mhd_LCLCEcntr_V2.m';
%
% iFcst = [2,3,6,7,8]
IFCST = [2,3,6,7,8];
%iFcst1 = 6;  % Forecast group to analyze, <10
%iFcst2 = 6;
itime1 = 1;  % 2011
itime2 = 2;  % 2012
irun1  = 1; 
irun2  = 7;

pthmat = pthmat1;

% F/cast info
load('hycom_tsis_expts.mat');  % EXPT array


ixx=0;
ymx = 0;
ifc=0;
for iFcst=IFCST
  ifc=ifc+1;
  for itime=itime1:itime2  % 2011 and 2012 time windows
    for irun=irun1:irun2
%      fmat1 = sprintf('%sMHD_LCLCE_nemo_persist_hycomfcst%2.2i-%2.2i%2.2i.mat',...
%                  pthmat,iFcst,itime,irun);
% Updated MHD calculation with LC/LCE combination that gives best MHD
      fmat1 = sprintf('%sMHD_LCLCE_nemo_hycomfcst%2.2i-%2.2i%2.2i-V2.mat',...
                  pthmat,iFcst,itime,irun);

      fprintf('Loading %s\n',fmat1);
      A=load(fmat1);  % MHD for f/cast and persistence
      DV=datevec(A.TM);
      ymax = max(max(A.MHD(:,1:2)));
%      ixx=ixx+1;
      MHD(ifc).TIME(itime).RUN(irun).mhd      = A.MHD(:,1);
      MHD(ifc).TIME(itime).RUN(irun).Name     = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
      MHD(ifc).TIME(itime).RUN(irun).Date_str = sprintf('%2.2i/%2.2i/%4.4i',DV(1,3:-1:1));
      MHD(ifc).TIME(itime).RUN(irun).mhdPrst  = A.MHD(:,2); % control NEMO true contour on day 1
%      MHD(ifc).hycom0= A.MHD(:,3); % HYCOM contour on day 1
      MHD(ifc).TIME(itime).RUN(irun).iFcst    = iFcst;  % F/cast group( h/cast initial fields)
      MHD(ifc).TIME(itime).RUN(irun).Hindcast_Number = iFcst;
      MHD(ifc).TIME(itime).RUN(irun).Time_period = itime;
      MHD(ifc).TIME(itime).RUN(irun).Run_number = irun;
      MHD(ifc).TIME(itime).RUN(irun).Ymax = ymax;
    end
  end
end
%Ntot = ixx;

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

% Persistence:
ipst = size(CLR,1)+1;
clrP = [0.4 0.4 0.4];
CLR(ipst,:) = clrP;

%MXD = A.MAXD;
TM=A.TM;
Td=TM-TM(1)+1;
Td=Td(:);
DV=datevec(TM);


% ----------------------
%
%  Plot time series of MHD
%  Select runs/f/casts to plot
% otherwise too many figures will open up
%  1 figure per 1 f-cast - 1run
% -----------------------
ifn = 0;
ixx = 0;
nfcst = 0;
ifc=0;
nRuns = length(MHD);
if f_mhd_tser == 1
  fprintf(' ===== Plotting MHD Time Series =====\n\n');
  ii=0;
  for iFcst=IFCST
    ii=ii+1;
    for itime=itime1:itime2  % 2011 and 2012 time windows
      ifn=ifn+1;
      for irun=irun1:irun1
        clr  = CLR(iFcst,:);
        mhd  = MHD(ii).TIME(itime).RUN(irun).mhd;
        mhdP = MHD(ii).TIME(itime).RUN(irun).mhdPrst;
        ymax = MHD(ii).TIME(itime).RUN(irun).Ymax;
								dstr = MHD(ii).TIME(itime).RUN(irun).Date_str;
								fcstname = MHD(ii).TIME(itime).RUN(irun).Name;
								stl = sprintf('MHD(km) HYCOM-NEMO Fcst# %i LCLCE cntrs, %s',iFcst,dstr);

        sub_plot_mhd_tser_v2(mhd,mhdP,clr,clrP,ifn,stl,ymax,irun);

								bottom_text(btx,'pwd',1,'Position',[0.09 0.02 0.4 0.04]);
      end
    end
  end

  fprintf('Time Series plotted ...\n');
  keyboard
end


% Weekly pools
WK=[1:7:92]; % weekly pools
nwk=length(WK)-1;

% Mean MHD by forecast/time groups
% For each forecast:
% Pool all runs by months/weeks together
ifc=0;
for iFcst = IFCST
  ifc=ifc+1;
  for itime=itime1:itime2
    POOL(ifc).Time(itime).pm1=[];  % month 1 pool
    POOL(ifc).Time(itime).pm2=[];
    POOL(ifc).Time(itime).pm3=[];

    for iwk=1:nwk
      POOL(ifc).Time(itime).WK(iwk).pw=[];
    end
  end
end

ifc=0;
for iFcst = IFCST
  ifc=ifc+1;
  for itime=itime1:itime2
    dmm = [];
    for irun=irun1:irun2
      mhd = MHD(ifc).TIME(itime).RUN(irun).mhd;
      mhdp = MHD(ifc).TIME(itime).RUN(irun).mhdPrst;

% Pool all runs for the same forecast group and time period
						pm1=POOL(ifc).Time(itime).pm1;
						pm2=POOL(ifc).Time(itime).pm2;
						pm3=POOL(ifc).Time(itime).pm3;

% Monthly pools
						pm1=[pm1;mhd(1:30)];
						pm2=[pm2;mhd(31:60)];
						pm3=[pm3;mhd(61:end)];

						POOL(ifc).Time(itime).pm1=pm1;
						POOL(ifc).Time(itime).pm2=pm2;
						POOL(ifc).Time(itime).pm3=pm3;

% Weekly pools
						for iwk=1:nwk
								id1=WK(iwk);
								id2=WK(iwk+1)-1;

								pw=POOL(ifc).Time(itime).WK(iwk).pw;
								pw=[pw;mhd(id1:id2)];
								POOL(ifc).Time(itime).WK(iwk).pw=pw;
						end
    end
  end
end

Nfgr = length(IFCST);
POOL(1).ylim=[0,160];
% ----------------------------
% Plot: all F/csts in 1 group for 3 30-day time periods
% coordinates of the Bars relative to X=# of the 30-day period, 1,2,3
nfg=40;
anls_nm='MHD(km) HYCOM LCLCE months';
sub_plotMHD_bars_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
bottom_text(btx,'pwd',1);

% Plot weekly bars
nfg=41;
anls_nm='MHD(km) HYCOM LCLCE weeks';
sub_plotMHD_bars_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,WK,anls_nm);
bottom_text(btx,'pwd',1);
   





