% Plot RMSE statistics
% T anom at 200 m
% Several OSSEs (IC from hindcast syntehic obs from 1km NEMO) and 
% Predictability 1a (initialized from interpolated NEMO)
% all 1a expts are pooled together
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear



pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';
btx = 'plot_RMSEbar_T200_OSSE_Prdct1a.m';


% Info about All hindcast experiments
load('hycom_tsis_expts.mat');

Ntime = 2;
IFCST=[ 3, 7, 8]; % forecast groups
Nfgr=length(IFCST); % # of forecasts groups 
Z0 = -200; 

%
% Combine RMSE fields
% 
ifc=0;
for ii=1:Nfgr
  iFcst=IFCST(ii);

  ifc=ifc+1;
  FCST = sub_fcst_info(iFcst);
  Nruns=1;
  nruns=Nruns;
  nmexp = FCST.Hindcast_Name;
  pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
  irun1 = FCST.run1;
  irun2 = FCST.run2;
  itime1= 1;
	itime2= 2;


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

      fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',iFcst,itime,irun);
	    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
  	  frmseout = sprintf('%sRMSE_ossefcst_dT%4.4i_%s.mat',pthout,abs(Z0),fcstname);

      clear RMSERR
      fprintf('Loading %s\n',frmseout);
      load(frmseout);

      ccn=ccn+1;
      dmm=RMSERR.ERR_squared;
      dmm=sqrt(nanmean(dmm)); % spatial average
      dmm=dmm(:);

      RMSE(ifc).Time(itime).Fcst_name=fcstname;
      RMSE(ifc).Time(itime).RMSE_mean(:,irun)=dmm;

    end % irun loop - 
  end
end;
% Weekly pools
WK=[1:7:92]; % weekly pools
nwk=length(WK)-1;

% Pool all runs for the same f/cast group
% Average by months:

for itime=1:2
  for ifc=1:Nfgr
    POOL(ifc).Time(itime).pm1=[];
    POOL(ifc).Time(itime).pm2=[];
    POOL(ifc).Time(itime).pm3=[];

    for iwk=1:nwk
      POOL(ifc).Time(itime).WK(iwk).pw=[];
    end

  end
end

for itime=1:Ntime
  for ifc=1:Nfgr
%
% Pool all runs - spatial mean RMSE
% groups by months
    dmm=RMSE(ifc).Time(itime).RMSE_mean;
    [a1,a2] = size(dmm);
    if isempty(dmm); continue; end;
    % Pool all runs for the same forecast group and time period
    pm1=POOL(ifc).Time(itime).pm1;
    pm2=POOL(ifc).Time(itime).pm2;
    pm3=POOL(ifc).Time(itime).pm3;

    A=dmm(1:30,:);
    [b1,b2]=size(A);
    A=reshape(A,[b1*b2,1]);
    pm1=[pm1;A];

    A=dmm(31:60,:);
    [b1,b2]=size(A);
    A=reshape(A,[b1*b2,1]);
    pm2=[pm2;A];

    A=dmm(61:a1,:);
    [b1,b2]=size(A);
    A=reshape(A,[b1*b2,1]);
    pm3=[pm3;A];

    POOL(ifc).Time(itime).pm1=pm1;
    POOL(ifc).Time(itime).pm2=pm2;
    POOL(ifc).Time(itime).pm3=pm3;

% Weekly pools
    for iwk=1:nwk
      id1=WK(iwk);
      id2=WK(iwk+1)-1;

      pw=POOL(ifc).Time(itime).WK(iwk).pw;
      A=dmm(id1:id2,:);
      [b1,b2]=size(A);
      A=reshape(A,[b1*b2,1]);
      pw=[pw;A];
      POOL(ifc).Time(itime).WK(iwk).pw=pw;
    end
  end
end

% Get Predict1a
POOL1a = sub_get_Prdct1aRMSE_T200;

ll=length(POOL);
for itime=1:2
  POOL(ll+1).Time(itime) = POOL1a.Time(itime);
end
IFCST(ll+1)=10;
EXPT(10).Name_short='Prdct1a P1-P7';
Nfgr = length(IFCST);

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

CLR(10,:)=CLR(2,:);

btx='plot_RMSEbar_T200_OSSE_Prdct1a.m';


% Plot: all F/csts in 1 group for 3 30-day time periods
% coordinates of the Bars relative to X=# of the 30-day period, 1,2,3
nfg=40;
anls_nm=sprintf('RMSE(m) Tanom 200m, months');
sub_plotMHD_bars_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
bottom_text(btx,'pwd',1);

% Plot weekly bars
nfg=41;
anls_nm=sprintf('RMSE(m) Tanom 200m, weeks');
sub_plotMHD_bars_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,WK,anls_nm);
bottom_text(btx,'pwd',1);







