% Predictability forecasts
%
% Analyze MHD for all forecasts
% 
% Before analysis:
% interpolate T to z0 vi     interp_fcst_hycom2z.m
% extract LC/LCE contour     extr_lc_tempPrdct_hycom.m
% compute MHD:               mhd_PrdctTz0cntr_nemo_hycom.m
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

T0 = 2.5;
Z0 = -200;

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

btx = 'anls_ALLfcstPrdct_mhd_Tz0LCLCEcntr.m';

IFCST=[10:16]; % Frecast groups to analyze
Nfgr=length(IFCST); % # of forecasts groups 

% Hindcast info
load('hycom_tsis_expts.mat');  % EXPT array


ixx=0;
ymx = 0;
for ii=1:Nfgr
  iFcst=IFCST(ii);

  FCST = sub_fcstPrdct_info(iFcst);
		Nhind = FCST.Nhind;  % initial cond from the  hindcast #
		hnd_name = FCST.Hindcast_Name;
		pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
		irun1 = FCST.run1;
  irun2 = FCST.run2;
  ntime = FCST.ntime; 
  Ntimes = ntime;

  for itime=1:ntime  % 2011 and 2012 time windows
    for irun=irun1:irun2
      nmexp   = sprintf('fcst%2.2i-%2.2i%2.2i',Nhind,itime,irun);
      fmat1 = sprintf('%sMHD_dT%2.2i_Z%3.3i_hycom_Prdct%s.mat',...
                 pthmat,round(T0*10),abs(Z0),nmexp);
      fprintf('Loading %s\n',fmat1);

      A=load(fmat1);  % MHD for f/cast and persistence
      DV=datevec(A.TM_fcst);
      ixx=ixx+1;
      MHD(ixx).mhd = A.MHD(:,1);
      MHD(ixx).Name = nmexp;
      MHD(ixx).Date_str = sprintf('%2.2i/%2.2i/%4.4i',DV(1,3:-1:1));
      MHD(ixx).nemo0 = A.MHD(:,2); % NEMO true contour on day 1 nemo persist
%      MHD(ixx).hycom0= A.MHD(:,2); % HYCOM contour on day 1 hycom persist
      MHD(ixx).iFcst=iFcst;    % Forecast group # (groupped by hindcast initial fields)
      MHD(ixx).Hindcast_Number=Nhind;
      MHD(ixx).Time_period=itime;
      MHD(ixx).Run_number=irun;
    end
  end
end
Ntot = ixx;

% Find corresponding hindcasts:
for ifc=1:Nfgr
  for ii=1:length(MHD);
    if MHD(ii).iFcst==IFCST(ifc); break; end;
  end
  HDCST(ifc)=MHD(ii).Hindcast_Number;
end

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
ipst = Ntot+1;
CLR(ipst,:) = [0.4 0.4 0.4];



%MXD = A.MAXD;
TM=A.TM_fcst;
Td=TM-TM(1)+1;
Td=Td(:);
DV=datevec(TM);


% Weekly pools
WK=[1:7:92]; % weekly pools
nwk=length(WK)-1;


% Mean MHD by forecast-hindcast/time groups
% Pool all mhd scores for same runs / months together
%
nll = length(MHD);
irun=0;
iFcst=[];
for ifc=1:Nfgr
  for itime=1:2
				POOL(ifc).Time(itime).pm1=[];
				POOL(ifc).Time(itime).pm2=[];
				POOL(ifc).Time(itime).pm3=[];
    for iwk=1:nwk
      POOL(ifc).Time(itime).WK(iwk).pw=[];
    end
  end
end

for ill=1:nll
  if isempty(iFcst); iFcst=MHD(ill).iFcst; end
  if MHD(ill).iFcst ~= iFcst;
    irun=0;
    iFcst=MHD(ill).iFcst;
  end
  irun=irun+1;
  ifc=find(IFCST==iFcst);
  itime=MHD(ill).Time_period;

% HYCOM mean MHD by 30 days
  dmm=MHD(ill).mhd;
%
% Pool all runs for the same forecast group and time period
  pm1=POOL(ifc).Time(itime).pm1;
  pm2=POOL(ifc).Time(itime).pm2;
  pm3=POOL(ifc).Time(itime).pm3;

  pm1=[pm1;dmm(1:30)];
  pm2=[pm2;dmm(31:60)];
  pm3=[pm3;dmm(61:end)];

  POOL(ifc).Time(itime).pm1=pm1;
  POOL(ifc).Time(itime).pm2=pm2;
  POOL(ifc).Time(itime).pm3=pm3;

% Weekly pools
  for iwk=1:nwk
    id1=WK(iwk);
    id2=WK(iwk+1)-1;

    pw=POOL(ifc).Time(itime).WK(iwk).pw;
    pw=[pw;dmm(id1:id2)];
    POOL(ifc).Time(itime).WK(iwk).pw=pw;
  end
end


%POOL(1).ylim=[0,160];
% Plot: all F/csts in 1 group for 3 30-day time periods
% coordinates of the Bars relative to X=# of the 30-day period, 1,2,3
nfg=40;
anls_nm=sprintf('MHD(km) \\DeltaT%3.1f z=%im, months',T0, Z0); 
sub_plotMHD_bars_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
bottom_text(btx,'pwd',1);

% Plot weekly bars
nfg=41;
anls_nm=sprintf('MHD(km) \\DeltaT%3.1f z=%im, weeks',T0, Z0);
sub_plotMHD_bars_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,WK,anls_nm);
bottom_text(btx,'pwd',1);


	






