% Tmeperature at 200 m
% Plot monthly mean maps of RMSE
% for 1 f/cast and 1 run
% for 2011 and 2012 
%
% Predictability experiments with
% experiments 1b - compare hycom perturbed f/casts 
% to the hycom control f/cast
% contolr run is done for 1a experiments (initialized from 
%  interpolated NEMO+GLORYS fields run for 100 days)
%
% Plot RMSE maps for forecast runs
% RMSE is averaged over the time periods: 1mo, 2 mo, 3 mo
% of the same forecast groups (i.e. initialized withe same hindcast)
%
% computed in calc_rmse_dTz0Prdct1b.m
%
%  Predictability experiments:
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

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

iFcst  = 13; % forecast #
% If need to average over all runs, irun1, ..., irun2, irun2> irun1
% For individual runs: irun2=irun1
% keep irun2=irun1
irun1  = 4; % irun =2, ..., 5 - for 
irun2  = irun1; 
itime1 = 1;  % May 2011
itime2 = 2;  % Jan 2012

Z0 = -200; % discard close to coastline points

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

FCST = sub_fcstPrdct_info(iFcst);

%
%Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
HH(isnan(HH))=100;
% Indices for subsampled GoM domain 
[LONs,LATs,HHs,INH,Iocn,xt1,xt2,yt1,yt2,its1,its2,jts1,jts2] = ...
   sub_subsample_HYCOM2GoM(HH,LON,LAT,Z0);


% 1 time period (2011/2012) at a time
Nruns=1; % # of forecast runs in each time window
nruns=Nruns;
iFcst = FCST.Nhind;  % initial cond from the  hindcast #
nmexp = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields

nfg=0;
for itime=itime1:itime2
%
% Combine RMSE fields
% 
  ccn=0;
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
    frmseout = sprintf('%sRMSE_Prdct1b_dT%4.4i_%s.mat',pthout,abs(Z0),fcstname);


    clear RMSERR
    fprintf('Loading %s\n',frmseout);
    load(frmseout);

    Iocn = RMSERR.Indx_gom;
    if ccn==0
      RMSE=Iocn*0; 
    end

    ccn=ccn+1;
    dmm=sqrt(RMSERR.ERR_squared);
    RMSE=RMSE+dmm;

  end % irun loop - 

  RMSE=RMSE/ccn;

  Time_name = sprintf('Start: %2.2i/%2.2i/%4.4i',DV(1,3),DV(1,2),DV(1,1));

  btx = 'plot_Tz0RMSE_Prdct1b_map.m';
% 
  for imo=1:3
% Average by months
    if imo==1
  				rmse1=nanmean(RMSE(:,1:31),2);  
    elseif imo==2
  				rmse1=nanmean(RMSE(:,32:61),2);  
    else
  				rmse1=nanmean(RMSE(:,62:90),2);  
    end

				nfg=nfg+1;
				stl = sprintf('Predct-1b Tz0Fcst# %i, run %2.2i: %s, Mo %i',iFcst,irun,Time_name,imo);
			
    c1=0;
    c2=2.5;
				sub_plot_rmse_v2(nfg,HHs,LONs,LATs,Iocn,rmse1,stl,c1,c2);
				bottom_text(btx,'pwd',1);
  end


end

