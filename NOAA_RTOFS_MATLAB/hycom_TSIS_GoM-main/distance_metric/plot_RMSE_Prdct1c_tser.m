% Plot RMSE statistics
% for forecast groups (i.e. initialized withe same hindcast)
%
% Predictability experiments with
% experiments 1 - 
% shift in atmospheric forcing
% control run and f/cast initialized at the same day
% atm focring is shifted by dt = +/- 1, 2 days
%
%
%  Predictability experiments:
%  For 2 time periods (May 2011 & Jan 2012)
%  NEMO+GLORYS interpolated into HYCOM
%      #10 - 01/05/2011 (Time 1) and 01/01/2012 (Time 2)
%      #11 - 08/05/2011 (Time 1) and 08/01/2012 (Time 2)
% etc
% Predictability experiments: iFcst 10, 11, ..., 16:
% Forecasts: use interpolated fields as initial fields
%            within each f/cast - perorm additional f/cast runs
%            initialized from day t0=7 days of HYCOM main f/cast
%            run1 = main f/cast, 
%            run2 = -2day shift used as IC at day t0
%            run3 = -1 day shift, etc.
%            run4 = +1 days
%            run5 = +2 days
% 
%


addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

% Plot 1 run for all f/casts at a time
%iFcst  = 13; % forecast #
irun1  = 5; % irun =2, ..., 5 - for 
irun2  = irun1;  % plot 1 irun (shifted IC) at a time
itime1 = 1;  % =1: May 2011
itime2 = 2;  % =2: Jan 2012
Ntimes = itime2-itime1+1;

f_prst = 1; % plot persistence RMSE 
Z0 = -200; % discard close to coastline points


DSHIFT = [0; -2; -1; 1; 2];

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';

% Info about All hindcast experiments
load('hycom_tsis_expts.mat');

%ts times:
FCST = struct;
cc=0;
dd=7;  % days shift in forecast start date
%ntime=1; % 2 time windows for forecasts
%irun1=1;
%irun2=Nruns;
imm=0;

%FCST = sub_fcstPrdct_info(iFcst);


%Ntime = 1;
IFCST = [10:16];
Nfgr=length(IFCST); % # of forecasts groups 

for ifc=1:Nfgr
  iFcst = IFCST(ifc);
  EXPT(iFcst).Name_short = sprintf('Prd1cFcst %2.2i',iFcst);
end



%
% Combine RMSE fields
% 
ifc=0;
for ii=1:Nfgr
  iFcst=IFCST(ii);

  ifc=ifc+1;
  FCST = sub_fcstPrdct_info(iFcst);
  Nruns=1;
		nruns=Nruns;
		iFcst = FCST.Nhind;  % initial cond from the  hindcast #
		nmexp = FCST.Hindcast_Name;
		pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
%		irun1 = FCST.run1;
%		irun2 = FCST.run2;
%  ntime = FCST.ntime;
%  ntime = 2;

  RMSE(ifc).Hnd_name=EXPT(iFcst).Name_short;
  RMSE(ifc).Fcst_nmb=iFcst;
  RMSE(ifc).Time(1).RMSE_mean=[];
  RMSE(ifc).Time(2).RMSE_mean=[];
  RMSE(ifc).Time(1).RMSEprs_mean = [];
  RMSE(ifc).Time(2).RMSEprs_mean = [];

		for itime=1:Ntimes
				ccn=0;
    RMSE(ifc).Time(itime).RMSE_mean=[];
				for irun=irun1:irun2    % forecast runs
						RUN = FCST.TIME0(itime).RUN(irun);
						pthd1 = RUN.pthbin;
						TM    = RUN.TM;
						YDAY  = RUN.jday;
						nrc   = length(TM);
						DV    = datevec(TM);

						fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',iFcst,itime,irun);
		%				fprintf(' Input data: %s\n',pthd1);
      fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
      frmseout = sprintf('%sRMSE%4.4im_Prdct1c_%s.mat',pthout,abs(Z0),fcstname);

						clear RMSERR
						fprintf('Loading %s\n',frmseout);
						load(frmseout);

						ccn=ccn+1;
						dmm=RMSERR.ERR_squared;
      dmm=sqrt(nanmean(dmm)); % spatial RMSE
      dmm=dmm(:);

      RMSE(ifc).Time(itime).Fcst_name=fcstname;
      RMSE(ifc).Time(itime).RMSE_mean(:,ccn)=dmm;
%
% Persistence
      pmm=RMSERR.ERRprst_squared;
      pmm=sqrt(nanmean(pmm)); % spatial RMSE persistence
      pmm=pmm(:);

      RMSE(ifc).Time(itime).RMSEprs_mean(:,ccn)=pmm; 
      
				end % irun loop - 

		end
end;
% Weekly pools
WK=[1:7:90]; % weekly pools
nwk=length(WK)-1;

% Pool all runs for the same f/cast group
% Average by months:

for itime=1:Ntimes
  for ifc=1:Nfgr
    POOL(ifc).Time(itime).pm1=[];
    POOL(ifc).Time(itime).pm2=[];
    POOL(ifc).Time(itime).pm3=[];

    POOL(ifc).Time(itime).pm1prs=[];
    POOL(ifc).Time(itime).pm2prs=[];
    POOL(ifc).Time(itime).pm3prs=[];

    for iwk=1:nwk
      POOL(ifc).Time(itime).WK(iwk).pw=[];
      POOL(ifc).Time(itime).WK(iwk).pw_prs=[];
    end

  end
end

for itime=1:Ntimes
  for ifc=1:Nfgr
%
% Pool all runs - spatial mean RMSE
% groups by months
				dmm=RMSE(ifc).Time(itime).RMSE_mean; % assuming 1 irun at a time
    if isempty(dmm); continue; end;
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

% Save time series:
    POOL(ifc).Time(itime).tser = dmm;
%
% Persistence:
    dmm=RMSE(ifc).Time(itime).RMSEprs_mean;
    if isempty(dmm); continue; end;
  % Pool all runs for the same forecast group and time period
    pm1=POOL(ifc).Time(itime).pm1prs;
    pm2=POOL(ifc).Time(itime).pm2prs;
    pm3=POOL(ifc).Time(itime).pm3prs;

    pm1=[pm1;dmm(1:30)];
    pm2=[pm2;dmm(31:60)];
    pm3=[pm3;dmm(61:end)];

    POOL(ifc).Time(itime).pm1prs=pm1;
    POOL(ifc).Time(itime).pm2prs=pm2;
    POOL(ifc).Time(itime).pm3prs=pm3;

  % Weekly pools
    for iwk=1:nwk
      id1=WK(iwk);
      id2=WK(iwk+1)-1;

      pw=POOL(ifc).Time(itime).WK(iwk).pw_prs;
      pw=[pw;dmm(id1:id2)];
      POOL(ifc).Time(itime).WK(iwk).pw_prs=pw;
    end
%
% Save time series:
    POOL(ifc).Time(itime).tser_prs = dmm;
  end
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


btx='plot_RMSE_Prdct1b_tser.m';

dsh = DSHIFT(irun);
if dsh>0
  srun = sprintf('Prd1c Run=%i (+%i day)',irun,abs(dsh));
else
  srun = sprintf('Prd1c Run=%i (-%i day)',irun,abs(dsh));
end

% Plot: all F/csts in 1 group for 3 30-day time periods
% coordinates of the Bars relative to X=# of the 30-day period, 1,2,3
nfg=10;
anls_nm=sprintf('%s RMSE(m) SSH mnth z>%im',srun, abs(Z0));
if f_prst==1
  sub_plotMHDprst_bars_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
else
  sub_plotMHD_bars_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
end
bottom_text(btx,'pwd',1);

% Plot weekly bars
nfg=11;
anls_nm=sprintf('%s RMSE(m) SSH weeks z>%im',srun,abs(Z0));
if f_prst==1
  sub_plotMHDprst_bars_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,WK,anls_nm);
else
  sub_plotMHD_bars_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,WK,anls_nm);
end
bottom_text(btx,'pwd',1);

% Plot t/series
nfg=12;
anls_nm=sprintf('%s RMSE(m) SSH days z>%im',srun,abs(Z0));
if f_prst==1
  sub_plotRMSEall_prst_tser(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm); % plot persistence as well
else
  sub_plotRMSEall_tser(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
end
bottom_text(btx,'pwd',1);








