% Load RMSE for PRedict. 1a experiments 
% T anomaly at 200 m depth
% and pool them together for 2011 and 2012
function POOL = sub_get_Prdct1aRMSE;

Z0 = -200;  % RMSE inside this depth 
f_prst = 1; % plot persistence RMSE 

irun1 = 1;
irun2 = 1; % there is only 1 run for each 1a experiment
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
%       fprintf(' Input data: %s\n',pthd1);

      fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
      frmseout = sprintf('%sRMSE_Prdct_dT%4.4i_%s.mat',pthout,abs(Z0),fcstname);

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
    POOL.Time(itime).pm1=[];
    POOL.Time(itime).pm2=[];
    POOL.Time(itime).pm3=[];

    for iwk=1:nwk
      POOL.Time(itime).WK(iwk).pw=[];
      POOL.Time(itime).WK(iwk).pw_prs=[];
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
    pm1=POOL.Time(itime).pm1;
    pm2=POOL.Time(itime).pm2;
    pm3=POOL.Time(itime).pm3;
%keyboard
    pm1=[pm1;dmm(1:30)];
    pm2=[pm2;dmm(31:60)];
    pm3=[pm3;dmm(61:91)];

    POOL.Time(itime).pm1=pm1;
    POOL.Time(itime).pm2=pm2;
    POOL.Time(itime).pm3=pm3;

    % Weekly pools
    for iwk=1:nwk
      id1=WK(iwk);
      id2=WK(iwk+1)-1;

      pw=POOL.Time(itime).WK(iwk).pw;
      pw=[pw;dmm(id1:id2)];
      POOL.Time(itime).WK(iwk).pw=pw;
    end

% Save time series:
    POOL.Time(itime).tser = dmm;
%
% Save time series:
    POOL.Time(itime).tser_prs = dmm;
  end
end

return

