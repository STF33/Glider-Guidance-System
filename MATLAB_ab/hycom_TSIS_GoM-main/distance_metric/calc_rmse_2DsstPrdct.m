% SST RMSE
% Layer 1 in HYCOM 3 m thick
% NEMO: Layers 1 = -1.05 m, Layer 3 (bottom interface) = -3.26198
%
% Predictability experiments:
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
%
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

%  Predictability forecasts (initialized from NEMO+GLORYS interpolation):
%                 #10 - NEMO+GLORYS at day 1
%           
iFcst = 10; % iFcst = hindcast # used for initial fields
f_mat = 2; % =1 save and overide existing mat; =2 - load saved and finish missing dates
rlrs_nemo = 3; % average 3 top layers for better agreement with HYCOM

if f_mat==0
  fprintf('Output is not saved !!! \n\n');
  fprintf('Check f_mat=%i\n',f_mat);
end

Z0 = -10; % discard close to coastline points

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';
pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';

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
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
fprintf('Loading NEMO grid %s\n',fgrd);
load(fgrd);
dZN = abs(diff(ZZN)); 



% 1 forecast group at a time
Nhnd = FCST.Nhind;  % initial cond from the  hindcast #
nmexp = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
irun1 = FCST.run1;
irun2 = FCST.run2;
ntime = FCST.ntime;
for itime=1:ntime
  for irun=irun1:irun2    % forecast runs
    RUN = FCST.TIME0(itime).RUN(irun);
    pthd1 = RUN.pthbin;
    TM    = RUN.TM;
    YDAY  = RUN.jday;
    nrc   = length(TM);
    DV    = datevec(TM);

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
    frmseout = sprintf('%sRMSEsst%4.4i_%s.mat',pthout,abs(Z0),fcstname);

    fprintf('\n\n %s %s\n',nmexp,frmseout);
				fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',Nhnd,itime,irun);
				fprintf(' Input data: %s\n',pthd1);

    clear RMSERR
    RMSERR.Name=EXPT(iFcst).Name;
    RMSERR.Name_short=EXPT(iFcst).Name_short;
    RMSERR.Indx_hycom=Iocn;
    RMSERR.HYCOM_dim=[m,n];
%    RMSERR.HYCOM_pnts=Iocn;

% Get persistence - day 1 from NEMO:
    dmean=0;
    iz0=1;
    dnmb=TM(1);
    ssm = 0;
    for iz0=1:rlrs_nemo
      sstn = sub_get_NEMO_TS(dnmb,'toce',iz0);
      ssm = ssm+dZN(iz0)*sstn;
    end
    sstn = ssm./sum(dZN(1:rlrs_nemo));
    Inm = find(INN==1);
    sstnM = nanmean(sstn(Inm));
    sstP = sstn-sstnM;  % demean
%keyboard
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
						iday= YDAY(ii);

% Day 1 - initial field - take from the hindcast:
						fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
						finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
      rlayer=1;
      fld = 'temp';
      [F,nn,mm,ll] = read_hycom(fina,finb,fld,'r_layer',rlayer);
      F(F>huge)=nan;
      sst=squeeze(F);
      I=find(INH==1);
      sstM=nanmean(sst(I));
      sst = sst-sstM; % demean


      dnmb=TM(ii);
%      sstn = sub_get_NEMO_TS(dnmb,'toce',iz0); 
						ssm = 0;
						for iz0=1:rlrs_nemo
								sstn = sub_get_NEMO_TS(dnmb,'toce',iz0);
								ssm = ssm+dZN(iz0)*sstn;
						end
						sstn = ssm./sum(dZN(1:rlrs_nemo));
      Inm = find(INN==1);
      sstnM = nanmean(sstn(Inm));
      sstn = sstn-sstnM;  % demean


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
        Err = [];
        ErrP= [];
        eH = sst(ih0);
        for itt=1:4
          in=Inemo(itt);
          jn=Jnemo(itt);
          Err(itt)=(eH-sstn(jn,in)).^2;
          ErrP(itt)=(sstP(jn,in)-sstn(jn,in)).^2; % persistence
        end

        cnn=cnn+1;
        RMSE(inx)=min(Err);
        RMSEP(inx)=min(ErrP);
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
        set(gca,'xlim',[1 600],...
                'ylim',[380 850]);
        keyboard
      end

      SSTM(ii) = sstM;
      SSTNM(ii) = sstnM;
      SSTM=SSTM(:);
      SSTNM=SSTNM(:);

      RMSE=RMSE(:);
      RMSERR.ERR_squared(:,ii)=RMSE;
      RMSERR.ERRprst_squared(:,ii)=RMSEP;  % persistence squared err
      RMSERR.MeanGoM_sst_hycom = SSTM;
      RMSERR.MeanGoM_sst_nemo = SSTNM;

      fprintf('1 record %6.3f min, RMSE=%6.4g, PerstRMSE=%6.4g, MnSST=%4.1f MnSST_nemo=%4.1f \n',...
              toc/60,rms,rmsp,sstM,sstnM);

      if mod(ii,15)==0 & f_mat==1  % f_mat=2 should change to f_mat=1
        fprintf('Saving %s\n',frmseout);
        save(frmseout,'RMSERR');
      end
    end

    if f_mat==1
      fprintf('Saving %s\n',frmseout);
      save(frmseout,'RMSERR');
    end

  end % irun loop - 
end

