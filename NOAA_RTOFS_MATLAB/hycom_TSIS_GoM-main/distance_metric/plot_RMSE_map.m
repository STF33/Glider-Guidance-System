% Plot RMSE maps for forecast runs
% RMSE is averaged over the time periods: 1mo, 2 mo, 3 mo
% of the same forecast groups (i.e. initialized withe same hindcast)
%
% computed in calc_rmse_2Dssh.m
% or parallel version calc_rmse_2Dssh_paral.m
%
% RMSE between NEMO and 
% HYCOM analysis SSH fields
% Forecasts: use HYCOM hindcasts for initial fields
%  Forecasts: decision tree:
%
%  Hindcast group (which is used to initialize f/cst): 3, 7 or 8 currently
% 2 - not finished
% H/cast #2 - Full 2D SSH  -- not finished
%        #3 - AVISO SSH tracks only
%        #6 - T/S GoM profiles - every 30th point on NEMO grid
%        #7 - AVISO + UGOS PIES (T/S profiles)
%        #8 - AVISO + extended PIES
%  Predictability experiments:
%        #10 - NEMO+GLORYS interpolated into HYCOM
%      |
%      V
%  TimePeriod during which the f/casts are initialized (May-June 2011, Jan-Feb 2012)
%     |
%     V
%  Forecasts (90 day f/csts shifted 7 days - 7 fcsts)
%
%
% Extract and save LC contours from HYCOM free runs / forecasts nd
% nemo simulations 
% specify individually which run need to extract
%
%  NEMO is extracted in extr_lc_hycom_nemo-V1.m 
%
% Compare LC contours from the HYCOM_TSIS hindcast
% assimilating NEMO 1/100 free running fields
% and NEMO LC contours
%


addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

iFcst = 3; % iFcst = hindcast # used for initial fields
itime0 = 2; % 2011/2012
Z0 = -10; % discard close to coastline points
%Z0 = -200; % discard shelf

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
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

FCST = sub_fcst_info(iFcst);


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

%XX=LON(jts1:jts2,its1:its2);
%YY=LAT(jts1:jts2,its1:its2);
%HHs=HH(jts1:jts2,its1:its2);
%Hnas=HHs;
%[msnsb]=size(HHs);




% 1 time period (2011/2012) at a time
Nruns=1; % # of forecast runs in each time window
nruns=Nruns;
Nhnd = FCST.Nhind;  % initial cond from the  hindcast #
nmexp = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
irun1 = FCST.run1;
irun2 = FCST.run2;
ntime = FCST.ntime;

for itime=itime0:itime0
%
% Combine RMSE fields
% 
  ccn=0;
  RMSE=Iocn*0; 
  for irun=irun1:irun2    % forecast runs
    RUN = FCST.TIME0(itime).RUN(irun);
    pthd1 = RUN.pthbin;
    TM    = RUN.TM;
    YDAY  = RUN.jday;
    nrc   = length(TM);
    DV    = datevec(TM);

    fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',Nhnd,itime,irun);
%    fprintf(' Input data: %s\n',pthd1);

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
%    frmseout = sprintf('%sRMSE_%s.mat',pthout,fcstname);    
    if Z0==-10
      frmseout = sprintf('%sRMSE_%s.mat',pthout,fcstname);
    else
      frmseout = sprintf('%sRMSE%4.4im_%s.mat',pthout,abs(Z0),fcstname);
    end

    clear RMSERR
    fprintf('Loading %s\n',frmseout);
    load(frmseout);

    ccn=ccn+1;
    dmm=sqrt(RMSERR.ERR_squared);
    RMSE=RMSE+dmm;

  end % irun loop - 

  RMSE=RMSE/ccn;

%  keyboard

  
  Hnd_name = EXPT(iFcst).Name_short;
  if itime==1
    Time_name='Apr/May 2011';
  else
    Time_name='Jan/Feb 2012';
  end

  btx = 'plot_RMSE_map.m';
% 
  for imo=1:3
% Average by months
    if imo==1
      rmse1=nanmean(RMSE(:,1:31),2);  
    elseif imo==2
      rmse1=nanmean(RMSE(:,32:61),2);  
    else
      rmse1=nanmean(RMSE(:,62:91),2);  
    end

    nfg=(itime-1)*3+imo;
    stl = sprintf('Init: %s, Fcsts# %i:  Z0=-%im %s, %i Mo',...
                   Hnd_name,iFcst,abs(Z0),Time_name,imo);
   
    sub_plot_rmse(nfg,HH,LON,LAT,Iocn,rmse1,stl);
    bottom_text(btx,'pwd',1);
  end


end

