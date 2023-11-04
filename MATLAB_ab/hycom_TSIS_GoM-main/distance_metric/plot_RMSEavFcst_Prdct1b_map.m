% Predictability experiments with
% experiments 1b - compare hycom perturbed f/casts 
% to the hycom control f/cast
% contolr run is done for 1a experiments (initialized from 
%  interpolated NEMO+GLORYS fields run for 100 days)
%
% Plot RMSE maps for given run (perturbation shift) averaged
% over all f/casts 10,..., 16 (inital dates)
%
% RMSE is averaged over the time periods: 1mo, 2 mo, 3 mo
% of the same forecast groups (i.e. initialized withe same hindcast)
%
% computed in calc_rmse_2DsshPrdct1b_hycom2hycom.m
%
% RMSE between NEMO and 
% HYCOM analysis SSH fields
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

% Plot averaged (over iFcst's) RMSE for each irun
irun   = 5; % irun =2, ..., 5 
% average over the forecasts:
ifcst1  = 10; % forecast #
ifcst2  = 16;
% for each time window:
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


%
% HYCOM:
rg=9806;  % convert pressure to depth, m
huge=1e20;

%Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
HH(isnan(HH))=100;


% GoM region HYCOM:
GOM=[366   489
   476   531
   583   560
   576   646
   508   827
   336   848
   204   829
    64   798
    19   746
    16   662
    12   578
    25   455
    71   382
   165   356
   281   400];

[XM,YM]=meshgrid([1:n],[1:m]);
INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
Iocn = find(INH==1 & HH<Z0);
Nocn=length(Iocn);
clear XM YM


% 
% Subsample to a smaller domain:
xnas1 = min(LON(Iocn));
xnas2 = max(LON(Iocn));
ynas1 = min(LAT(Iocn));
ynas2 = max(LAT(Iocn));
[J,I]=find(LON<xnas1);
its1=max(I);
[J,I]=find(LON>xnas2);
its2=min(I);
[J,I]=find(LAT<ynas1);
jts1=max(J);
[J,I]=find(LAT>ynas2);
jts2=min(J);

xt1=LON(jts1,its1);
yt1=LAT(jts1,its1);
xt2=LON(jts2,its2);
yt2=LAT(jts2,its2);



% 1 time period (2011/2012) at a time
%Nruns=1; % # of forecast runs in each time window
%nruns=Nruns;
%iFcst = FCST.Nhind;  % initial cond from the  hindcast #
%nmexp = FCST.Hindcast_Name;
%pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
%irun1 = FCST.run1;
%irun2 = FCST.run2;
%ntime = FCST.ntime;

nfg=0;
for itime=itime1:itime2
%
% Combine RMSE fields
% 
  ccn=0;
  RMSE=Iocn*0; 
  for iFcst=ifcst1:ifcst2    % forecast runs
    FCST = sub_fcstPrdct_info(iFcst);

    RUN = FCST.TIME0(itime).RUN(irun);
    pthd1 = RUN.pthbin;
    TM    = RUN.TM;
    YDAY  = RUN.jday;
    nrc   = length(TM);
    DV    = datevec(TM);

				fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',iFcst,itime,irun);
%				fprintf(' Input data: %s\n',pthd1);

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
    frmseout = sprintf('%sRMSE%4.4im_Prdct1b_%s.mat',pthout,abs(Z0),fcstname);


    clear RMSERR
    fprintf('Loading %s\n',frmseout);
    load(frmseout);

    ccn=ccn+1;
    dmm=sqrt(RMSERR.ERR_squared);
    RMSE=RMSE+dmm;

  end % irun loop - 

  RMSE=RMSE/ccn;

  Time_name = sprintf('Start: %2.2i/%2.2i/%4.4i',DV(1,3),DV(1,2),DV(1,1));

  btx = 'plot_RMSEavFcst_Prdct1b_map.m';
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
				stl = sprintf('Predct1b avrg:Fcst# %i-%i, run %2.2i: %s, Mo %i',ifcst1,ifcst2,irun,Time_name,imo);
			
				sub_plot_rmse(nfg,HH,LON,LAT,Iocn,rmse1,stl);
				bottom_text(btx,'pwd',1);
  end


end

