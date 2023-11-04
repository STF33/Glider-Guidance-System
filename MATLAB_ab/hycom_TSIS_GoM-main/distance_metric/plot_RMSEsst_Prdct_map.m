% Predictability experiments with
% NEMO_GLORYS interpolated as initial state
%
% Plot RMSE maps for SST (1st HYCOM layer)
% forecast runs
% RMSE is averaged over the time periods: 1mo, 2 mo, 3 mo
% of the same forecast groups (i.e. initialized withe same hindcast)
%
% computed in calc_rmse_2DsstPrdct.m
%
% RMSE between NEMO and 
% HYCOM analysis SSH fields
%
%  Predictability experiments:
%  NEMO+GLORYS interpolated into HYCOM
%      #10 - 01/05/2011 (Time 1) and 01/01/2012 (Time 2)
%      #11 - 08/05/2011 (Time 1) and 08/01/2012 (Time 2)
% etc
%      |
%      V
%  TimePeriod during which the f/casts are initialized (May-June 2011, Jan-Feb 2012)
%     |
%     V
%  Forecasts (100 day f/csts shifted 7 days - 7 fcsts)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

iFcst = 10; % iFcst = hindcast # used for initial fields

Z0 = -10; % discard close to coastline points

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

FCST = sub_fcstPrdct_info(iFcst);


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

for itime=1:ntime
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
%				fprintf(' Input data: %s\n',pthd1);

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
    frmseout = sprintf('%sRMSEsst%4.4i_%s.mat',pthout,abs(Z0),fcstname);    

    clear RMSERR
    fprintf('Loading %s\n',frmseout);
    load(frmseout);

    ccn=ccn+1;
    dmm=sqrt(RMSERR.ERR_squared);
    RMSE=RMSE+dmm;

    Thycom = RMSERR.MeanGoM_sst_hycom;
    Tnemo  = RMSERR.MeanGoM_sst_nemo;
    daa = Thycom-Tnemo;
    dT(:,irun) = daa(:);
  end % irun loop - 

  RMSE=RMSE/ccn;

%  keyboard

  
  Hnd_name = EXPT(iFcst).Name_short;
  Time_name = sprintf('Start: %2.2i/%2.2i/%4.4i',DV(1,3),DV(1,2),DV(1,1));

  btx = 'plot_RMSEsst_Prdct_map.m';
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
				stl = sprintf('Predctb RMSE SST Fcst# %i, run %2.2i: %s, Mo %i',iFcst,irun,Time_name,imo);
		
    c1=0;
    c2=4;	
    sub_plot_rmseSST(nfg,HH,LON,LAT,Iocn,rmse1,stl,c1,c2);
				bottom_text(btx,'pwd',1);
  end

%
% Plot GoM mean SST:
  ifg=nfg+10;
  figure(ifg); clf;
  set(gcf,'Position',[249 450 1036 805]);
  axes('Position',[0.09 0.55 0.86 0.35]);
  hold on;
  plot(dT,'-','Linewidth',2,'Color',[0 0.5 0.8]);
  stl = sprintf('SST_{hycom}-SST_{nemo}, Fcst# %i, run %2.2i: %s',iFcst,irun,Time_name);
  title(stl,'Fontsize',12);
  set(gca,'tickdir','out',...
          'xlim',[0 length(dT)],...
          'xtick',[0:10:150],...
          'xgrid','on',...
          'ygrid','on',...
          'Fontsize',14);
  xlabel('Days');
  bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.5 0.05]);



end

