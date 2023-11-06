% Plot RMSE statistics
% for forecast groups (i.e. initialized withe same hindcast)
%  Predictability 1a - compare HYCOM to NEMO
% HYCOM initialized from interpolated NEMO and ran for ~100 days
%
% computed in calc_rmse_2Dssh.m
% or parallel version calc_rmse_2Dssh_paral.m
%
% RMSE between NEMO and 
% HYCOM analysis SSH fields
% Forecasts: use HYCOM hindcasts for initial fields
%  Forecasts: decision tree:
%
%  Predictability experiments:
%  For 2 time periods (May 2011 & Jan 2012)
%  initalize forecasts and run for 100 days
%  shift f/casts by 7 days
%  initial fields - interpolated NEMO+GLORYS
% 
%        #10 - NEMO+GLORYS interpolated into HYCOM
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
startup;

close all
clear

%Z0=-10;   % ~whole GoM
Z0 = -200;  % RMSE inside this depth 
f_prst = 1; % plot persistence RMSE 

irun1 = 1; 
irun2 = 1; % only 1 run for each experiment
itime1 = 1; 
itime2 = 2;

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
%pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';

% Info about All hindcast experiments
load('hycom_tsis_expts.mat');

%Ntime = 1;
IFCST = [10:16];
Nfgr=length(IFCST); % # of forecasts groups 

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
	Nhnd = FCST.Nhind;  % initial cond from the  hindcast #
	nmexp = FCST.Hindcast_Name;
	pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields


  RMSE(ifc).Hnd_name=EXPT(iFcst).Name_short;
  RMSE(ifc).Fcst_nmb=iFcst;
  RMSE(ifc).Time(1).RMSE_mean=[];
  RMSE(ifc).Time(2).RMSE_mean=[];
  RMSE(ifc).Time(1).RMSEprs_mean = [];
  RMSE(ifc).Time(2).RMSEprs_mean = [];

	for itime=itime1:itime2
		ccn=0;
    RMSE(ifc).Time(itime).RMSE_mean=[];
		for irun=irun1:irun2    % forecast runs
			RUN = FCST.TIME0(itime).RUN(irun);
			pthd1 = RUN.pthbin;
			TM    = RUN.TM;
			YDAY  = RUN.jday;
			nrc   = length(TM);
			DV    = datevec(TM);

			fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',Nhnd,itime,irun);
%				fprintf(' Input data: %s\n',pthd1);

			fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
			if abs(Z0)>10
				frmseout = sprintf('%sRMSE%4.4im_%s.mat',pthout,abs(Z0),fcstname);
			else
				frmseout = sprintf('%sRMSE_%s.mat',pthout,fcstname);
			end


			clear RMSERR
			fprintf('Loading %s\n',frmseout);
			load(frmseout);

			ccn=ccn+1;
			dmm=RMSERR.ERR_squared;
      dmm=sqrt(nanmean(dmm)); % spatial RMSE
      dmm=dmm(:);

      RMSE(ifc).Time(itime).Fcst_name=fcstname;
      RMSE(ifc).Time(itime).RMSE_mean(:,irun)=dmm;
%
% Persistence
      pmm=RMSERR.ERRprst_squared;
      pmm=sqrt(nanmean(pmm)); % spatial RMSE persistence
      pmm=pmm(:);

      RMSE(ifc).Time(itime).RMSEprs_mean(:,irun)=pmm; 
      
		end % irun loop - 
	end
end;
% Weekly pools
WK=[1:7:92]; % weekly pools
nwk=length(WK)-1;

% Pool all runs for the same f/cast group
% Average by months:

for itime=itime1:itime2
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

for itime=itime1:itime2
  for ifc=1:Nfgr
%
% Pool all runs - spatial mean RMSE
% groups by months
		dmm=RMSE(ifc).Time(itime).RMSE_mean;
    if isempty(dmm); continue; end;
		% Pool all runs for the same forecast group and time period
		pm1=POOL(ifc).Time(itime).pm1;
		pm2=POOL(ifc).Time(itime).pm2;
		pm3=POOL(ifc).Time(itime).pm3;
%keyboard
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


btx='plot_RMSE_Predict_tser.m';

% Plot: all F/csts in 1 group for 3 30-day time periods
% coordinates of the Bars relative to X=# of the 30-day period, 1,2,3
nfg=10;
anls_nm=sprintf('RMSE(m) SSH mnth z>%im',abs(Z0));
if f_prst==1
  sub_plotMHDprst_bars_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
else
  sub_plotMHD_bars_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
end
bottom_text(btx,'pwd',1);

% Plot weekly bars
nfg=11;
anls_nm=sprintf('RMSE(m) SSH weeks z>%im',abs(Z0));
if f_prst==1
  sub_plotMHDprst_bars_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,WK,anls_nm);
else
  sub_plotMHD_bars_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,WK,anls_nm);
end
bottom_text(btx,'pwd',1);

% Plot t/series
nfg=12;
anls_nm=sprintf('RMSE(m) SSH days z>%im',abs(Z0));
if f_prst==1
  sub_plotRMSEall_prst_tser(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm); % plot persistence as well
else
  sub_plotRMSEall_tser(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
end
bottom_text(btx,'pwd',1);








