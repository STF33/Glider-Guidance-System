% Predictability experiments:  hycom-to-hycom comparison
%
% experiments 1c - perturbation in atmospheric forcing
%                  use perturbed runs from 1b 
%                  but compare with the control run
%                  that starts from the same day
%                  same initial conditions but different atm forcing
%    i.e. dt=-2 days: perturbed run starts at t0-2 but uses
%    atm forcing from t0, whereas control run
%    starts at t0-2 and uses atm forcing from t0-2 (correct forcing)
%    
% hycom control f/cast
% contolr run is done for 1a experiments (initialized from 
%  interpolated NEMO+GLORYS fields run for 100 days)
%
% f/casts initialized from interpolated NEMO+GLORYS
%
% Calculate RMSE between NEMO and 
% HYCOM analysis SSH fields
% Predictability experiments: iFcst 10, 11, ..., 16:
% Forecasts: use interpolated fields as initial fields
%            within each f/cast - perorm additional f/cast runs
%            initialized from day t0=7 days of HYCOM main f/cast
%            run1 = main f/cast, 
%            run2 = -2day shift used as IC at day t0
%            run3 = -1 day shift, etc.
%			         run4 = +1 days
%            run5 = +2 days
%
%


addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

%  Predictability forecasts (initialized from NEMO+GLORYS interpolation):
%                 #10 - NEMO+GLORYS at day 1
%           
iFcst = 16; % iFcst = hindcast # used for initial fields
irun1=2;   % irun1 - control f/cast, 2,.., 5 - shifted +/- 2, 1 perturbation f/casts
irun2=5;
itime1=1;   % =1 May 2011, 2 = Jan 2012
itime2=2;   % 

f_mat = 1; % =1 save and overide existing mat; =2 - load saved and finish missing dates

if f_mat==0
  fprintf('\nOutput is not saved !!! \n\n');
  fprintf('Check f_mat=%i\n',f_mat);
end

%
% Specify depth limit for RMSE calculation
%Z0 = -10; % discard close to coastline points
Z0 = -200; % only deep GoM

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';



% Info about All hindcast experiments
load('hycom_tsis_expts.mat');


%ts times:
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



%keyboard

% 


% 1 forecast group at a time
Nhnd = FCST.Nhind;  % initial cond from the control run
nmexp = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
for itime=itime1:itime2
  for irun=irun1:irun2    % forecast runs
    RUN0  = FCST.TIME0(itime).RUN(1); % control run
    RUN   = FCST.TIME0(itime).RUN(irun);
    pthd1 = RUN.pthbin;
    TM    = RUN.TM;
    YDAY  = RUN.jday;
    nrc   = length(TM);
    day_shift = RUN.dayT0_true-RUN.prdct_time0;
    TM    = TM+day_shift;
    DV    = datevec(TM);

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
    frmseout = sprintf('%sRMSE%4.4im_Prdct1c_%s.mat',pthout,abs(Z0),fcstname);

    fprintf('\n\n %s %s\n',nmexp,frmseout);
				fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',Nhnd,itime,irun);
				fprintf(' Input data: %s\n',pthd1);
    fprintf(' Output saved: %s\n',frmseout);

    clear RMSERR
    RMSERR.Name=EXPT(iFcst).Name;
    RMSERR.Name_short=EXPT(iFcst).Name_short;
    RMSERR.Indx_hycom=Iocn;
    RMSERR.HYCOM_dim=[m,n];
%    RMSERR.HYCOM_pnts=Iocn;

% Get persistence 
% Persistence - control run from HYCOM control run
    dmean=1;
% Time t0:
    dnmb0 = FCST.TIME0(itime).RUN(irun).dayT0_true;
    dv0 = datevec(dnmb0);
    pthd_cntrl = FCST.TIME0(itime).RUN(1).pthbin;
    yr = dv0(1);
    jday0 = dnmb0-datenum(yr,1,1)+1;
    if dnmb0~=TM(1),
      error('Mismatched day 1 for control run and perturb exprt');
    end
    
    fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd_cntrl,yr,jday0);
    finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd_cntrl,yr,jday0);
    sshP = sub_getSSH_hycom(fina,finb,INH,dmean);
    
% Continue from the last saved record
% skip finished mat files
    last_rec=0;
    if f_mat==2
      fprintf('Continue from the last unfinished mat file\n');
      if exist(frmseout,'file');
        fprintf('Loading %s\n',frmseout);
        load(frmseout);
        [a1,a2]=size(RMSERR.ERR_squared);
        last_rec=a2;
%        keyboard
      end
    end
        

%				cntr=0;
				for ii=1:nrc
      if f_mat==2
        if ii<=last_rec 
          fprintf('%s Record exist skipping %i/%i/%i\n',fcstname,DV(ii,1:3));
          continue
        else
          fprintf('Continue from ii=%i %s %i/%i/%i\n\n',ii,fcstname,DV(ii,1:3));
          f_mat=1;
          last_rec=0;
        end
      end
%      keyboard
						tic;

						yr  = DV(ii,1);
						mo  = DV(ii,2);
						dm  = DV(ii,3);
						dnmb= TM(ii);
						iday= YDAY(ii); % year days in the perturbed simulations, shifted days
      jday_shift = iday+day_shift; % true year day

% Day 1 - initial field - take from the control run:
      fprintf('Reading f/cast \n');
      if dnmb==dnmb0
        fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd_cntrl,yr,jday_shift);
        finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd_cntrl,yr,jday_shift);
      else
								fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
								finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
      end
      ssh = sub_getSSH_hycom(fina,finb,INH,dmean);

% Control run:
      fprintf('Reading control run \n');
      fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd_cntrl,yr,jday_shift);
      finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd_cntrl,yr,jday_shift);
      sshn = sub_getSSH_hycom(fina,finb,INH,dmean);

      RMSE=[];
      RMSE  = (ssh(Iocn)-sshn(Iocn)).^2;  % f/cast - control run
      RMSEP = (sshP(Iocn)-sshn(Iocn)).^2;  % Persistence - control run

%keyboard
      rms=sqrt(nansum(RMSE)/Nocn);
      rmsp=sqrt(nansum(RMSEP)/Nocn);

% Spatial map of RMSE;
      f_chck=0;
      if f_chck==1
        ERR=zeros(m,n)*nan;
        ERR(Iocn)=RMSE;
        figure(11); clf;
        pcolor(ERR); shading flat;
      end

      RMSE=RMSE(:);
      RMSERR.ERR_squared(:,ii)=RMSE;
      RMSERR.ERRprst_squared(:,ii)=RMSEP;  % persistence squared err

      fprintf('1 record %6.3f min, RMSE=%8.4g, PerstRMSE=%8.4g \n',toc/60,rms,rmsp);

      if mod(ii,20)==0 & f_mat==1  % f_mat=2 should change to f_mat=1
        fprintf('Saving %s\n',frmseout);
        save(frmseout,'RMSERR');
      end
%      keyboard
    end

    if f_mat==1
      fprintf('Saving %s\n',frmseout);
      save(frmseout,'RMSERR');
    end

  end % irun loop - 
end

