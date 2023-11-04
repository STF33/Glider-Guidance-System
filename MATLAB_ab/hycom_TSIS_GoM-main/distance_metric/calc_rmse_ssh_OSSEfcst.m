% Calculate RMSE between NEMO and 
% HYCOM SSH fields from OSSE forecasts 
% Forecasts: use HYCOM hindcasts for initial fields
%
%  Hindcast group (which is used to initialize f/cst): 3, 7 or 8 currently
% 2 - not finished
% H/cast #2 - Full 2D SSH
%        #3 - AVISO SSH tracks only
%        #6 - T/S profiles GoM, every 30th pnt NEMO
%        #7 - AVISO + UGOS PIES (T/S profiles)
%        #8 - AVISO + extended PIES
%      |
%      V
%  TimePeriod during which the f/casts are initialized (May-June 2011, Jan-Feb 2012)
%     |
%     V
%  Forecasts (90 day f/csts shifted 7 days - 7 fcsts)
%
% Predictability experiments: iFcst 10, 11, ..., 16:
%  use calc_rmse_2DsshPrdct.m 
%


addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

% OSSEs with NEMO subsampled for observed locations:
% iFcst =   H/cast #2 - Full 2D SSH   - only 2011
%                 #3 - AVISO SSH tracks only
%                 #6 - T/S GoM at every 30th pnt on NEMO grid - only 2011
%                 #7 - AVISO + UGOS PIES (T/S profiles)
%                 #8 - AVISO + extended PIES
%  Predictability forecasts (initialized from NEMO+GLORYS interpolation) 1a and 1a:
%           use calc_rmse_2DsshPrdct*
%
% 
iFcst = 8; % iFcst = hindcast # used for initial fields
f_mat = 2; % save mat; =2 - load saved and finish missing dates

%Z0 = -10; % discard close to coastline points
Z0 = -200;

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';

% NEMO 2 HYCOM indices:
%fmatout = sprintf('%snemo2hycom_indx.mat',pthmat);
%load(fmatout);


% Info about All hindcast experiments
load('hycom_tsis_expts.mat');

% Load NEMO 2 hycom indices:
fmatindx = sprintf('%snemo2hycom_indx.mat',pthmat);
fprintf('Loading %s\n',fmatindx);
load(fmatindx);

%ts times:
FCST = struct;
cc=0;
dd=7;  % days shift in forecast start date
Nruns=7; % # of forecast runs in each time window
ntime=2; % 2 time windows for forecasts
irun1=1;
irun2=Nruns;
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


if m~=INDX.HYCOM_size(1) | n~=INDX.HYCOM_size(2)
  fprintf('Saved domain: %i %i, this domain: %i %i\n',...
   INDX.HYCOM_size(1:2), m, n)
  error('Saved NEMO indices are for different HYCOM domain');
end


IJhycom = INDX.HYCOM_ocean;  % saved HYCOM pnts with NEMO indices


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
[LONN,LATN,INN] = sub_getNEMO_info;



% 1 forecast group at a time
nruns=Nruns;
Nhnd = FCST.Nhind;  % initial cond from the  hindcast #
nmexp = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
irun1 = FCST.run1;
irun2 = FCST.run2;
ntime = FCST.ntime;
%for itime=1:ntime
%  for irun=irun1:irun2    % forecast runs
for itime=1:1
  for irun=5:5
    RUN = FCST.TIME0(itime).RUN(irun);
    pthd1 = RUN.pthbin;
    TM    = RUN.TM;
    YDAY  = RUN.jday;
    nrc   = length(TM);
    DV    = datevec(TM);

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
    frmseout = sprintf('%sRMSE_%s.mat',pthout,fcstname);

    fprintf('\n\n %s %s\n',nmexp,frmseout);
    fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',Nhnd,itime,irun);
    fprintf(' Input data: %s\n',pthd1);

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
    frmseout = sprintf('%sRMSE_%s.mat',pthout,fcstname);   
    if abs(Z0)>10
      frmseout = sprintf('%sRMSE%4.4im_%s.mat',pthout,round(abs(Z0)),fcstname); 
    end

    clear RMSERR
    RMSERR.Name=EXPT(iFcst).Name;
    RMSERR.Name_short=EXPT(iFcst).Name_short;
    RMSERR.Indx_hycom=Iocn;
    RMSERR.HYCOM_dim=[m,n];
%    RMSERR.HYCOM_pnts=Iocn;


    dmean=1;
% Get persistence - day 1 from NEMO:
%    sshP = sub_getSSH_nemo(DV(1,:),LONN,LATN,INN,dmean);
% Persistence - 1st day of f.cast:
    ii = 1;
    yr  = DV(ii,1);
    mo  = DV(ii,2);
    dm  = DV(ii,3);
    dnmb= TM(ii);
    iday= YDAY(ii);
    fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthHcst,yr,iday);
    finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthHcst,yr,iday);
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
        

%    cntr=0;
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
      iday= YDAY(ii);

% Day 1 - initial field - take from the hindcast:
      if ii==1
        fprintf('Initial fields from hindcast\n');
        fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthHcst,yr,iday);
        finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthHcst,yr,iday);
      else
        fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
        finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
      end
      dmean=1;
      ssh = sub_getSSH_hycom(fina,finb,INH,dmean);

      sshn = sub_getSSH_nemo(DV(ii,:),LONN,LATN,INN,dmean);

      RMSE=[];
      cnn=0; 
      for inx=1:Nocn
        if mod(inx,10000)==0,
          fprintf(' Calculating RMSE: %4.1f%% done ...\n',inx/Nocn*100);
        end

        ih0=Iocn(inx);
        II=find(IJhycom==ih0);
        if isempty(II); 
          RMSE(inx)=nan;
          continue; 
        end; % pnt outside saved indices

        Inemo=INDX.I_NEMO(II,:);
        Jnemo=INDX.J_NEMO(II,:);

        f_chck=0;
        if f_chck==1
          figure(10); clf;
          hold on;
          plot(LON(ih0),LAT(ih0),'o');

          for itt=1:4
            if itt>1
            plot(LONN(Jnemo(itt),Inemo(itt)),LATN(Jnemo(itt),Inemo(itt)),'k*');
            else
            plot(LONN(Jnemo(itt),Inemo(itt)),LATN(Jnemo(itt),Inemo(itt)),'r*');
            end
          end
        end

%
% Take the smallest RMSD
% out of the surrounding closest NEMO points:
        Err  = [];
        ErrP = [];
        eH   = ssh(ih0);
        eHp  = sshP(ih0);
        for itt=1:4
          in=Inemo(itt);
          jn=Jnemo(itt);
          Err(itt)=(eH-sshn(jn,in)).^2;
%          ErrP(itt)=(sshP(jn,in)-sshn(jn,in)).^2; % persistence
          ErrP(itt)=(eHp-sshn(jn,in)).^2; % persistence
        end

        cnn=cnn+1;
        RMSE(inx)=min(Err);
        RMSEP(inx)=min(ErrP);
%        Imin(inx)=find(Err==min(Err),1);

      end;  % ocean points
%keyboard
      rms=sqrt(nansum(RMSE)/cnn);
      rmsp=sqrt(nansum(RMSEP)/cnn);

% Spatial map of RMSE;
      f_chck=0;
      if f_chck==1
        ERR=zeros(m,n)*nan;
        ERR(Iocn)=sqrt(RMSE);
        figure(11); clf;
        pcolor(ERR); shading flat;
      end

      RMSE=RMSE(:);
      RMSERR.ERR_squared(:,ii)=RMSE;
      RMSERR.ERRprst_squared(:,ii)=RMSEP;  % persistence squared err

%      fprintf('1 record %6.3f min, RMSE=%8.4f \n',toc/60,rms);
      fprintf('1 record %6.3f min, RMSE=%8.4g, PerstRMSE=%8.4g \n',toc/60,rms,rmsp);

      if mod(ii,15)==0 & f_mat==1
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

