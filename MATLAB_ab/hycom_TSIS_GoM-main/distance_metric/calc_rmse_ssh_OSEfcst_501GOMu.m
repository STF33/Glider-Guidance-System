% Calculate RMSE between OSE (PIES or noPIES)
% in the OSE forecasts (3 months f/casts)
% started every month for PIES and noPIES
% May 1 2009 - Dec 1 2010
% and 50.1GOMl reanalysis
% The f/casts are 
% Initialized from OSE hindcasts: PIES/noPIES
%
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

YR1 = 2010;
YR2 = YR1;
esimH = 'PIES';
%esimH = 'noPIES';
z0 = -200; 

f_mat = 1; % =1 save and overide existing mat; =2 - load saved and finish missing dates

if f_mat==0
  fprintf('\nOutput is not saved !!! \n\n');
  fprintf('Check f_mat=%i\n',f_mat);
end

%
% Specify depth limit for RMSE calculation
Z0 = -200; % only deep GoM

pthtopo      = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat       = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout       = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst     = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';
pthmat_sshR  = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthtopo_rnls = '/nexsan/archive/GOMu0.04_501/topo/';


% HYCOM:
rg=9806;  % convert pressure to depth, m
huge=1e20;

%Read 0.03 IAS HYCOM topography:
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

% IAS - GOMu indices
% see find_gomu2ias_indx.m
fmatindx = sprintf('%sgomu2ias_indx.mat',pthmat);
fprintf('Loading %s\n',fmatindx);
load(fmatindx);
%IX_gom = INDX.IAS_HYCOM_Iocn;  % linear IAS indices for saved IAS-GOMu mapping
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


%
% Topo, grid for 501 GOMu:
pthtopoR = '/nexsan/archive/GOMu0.04_501/topo/';
ftopoR = sprintf('%sdepth_GOMu0.04_03i.nc',pthtopoR);
HHR  = -1*squeeze(nc_varget(ftopoR,'depth'));
HHR(isnan(HHR))=100;
latR = nc_varget(ftopoR,'Latitude');
lonR = nc_varget(ftopoR,'Longitude');
[mR,nR]=size(HHR);
HHR(isnan(HHR))=100;

I = find(lonR>180.);
lonR(I) = lonR(I)-360.;

[LONR,LATR] = meshgrid(lonR,latR);


GOMR = sub_find_501GOMindx(GOM,LON,LAT,lonR,latR);
[XM,YM] = meshgrid([1:nR],[1:mR]);
INR = inpolygon(XM,YM,GOMR(:,1),GOMR(:,2));
Irnl = find(INR==1 & HHR<Z0);
clear XM YM

SSHM = [];
for YR=YR1:YR2
  ys = YR;
  iyr = YR;

  im1=1;
  im2=12;
  if iyr==2009, im1=5; end;

%  for imo=im1:im2
  for imo=im1:11
% Forecast : 3 months
% Dates for forecast:

% To fix error - skip some months:
    fprintf('Year %i, month %i\n',YR,imo);
    if iyr == 2009 & imo<10, continue; end;
    if iyr == 2010 & imo<10, continue; end;
    

    jd1 = datenum(iyr,1,1);
    id11 = datenum(iyr,imo,1)-jd1+1; % start Yr. day
    dnmb1 = datenum(iyr,imo,1);
 
    dnmb2=datenum(ys,imo,1)+100;
    dv=datevec(dnmb2);
    ye=dv(1);
    ime=dv(2);
    dnmb2=datenum(ye,ime,1)-1;
    jd2=datenum(ye,1,1);
    dv2=datevec(dnmb2);
    ye=dv2(1);
    id22=dnmb2-jd2+1;
    ndays = dnmb2-dnmb1+1;

    cc=0;
    dnmb=dnmb1-1;
    YRPLT=[];
    for idd=1:ndays
% Jan 1 is missing and started from Jan 2
      dnmb=dnmb+1;
      if imo==1 & dnmb==dnmb1,
        continue;
      end
      cc=cc+1;
      DV=datevec(dnmb);
      jd1=datenum(DV(1),1,1);
      yday=dnmb-jd1+1;
      YRPLT(cc,1)=DV(1);
      YRPLT(cc,2)=yday;
      YRPLT(cc,3)=dnmb;
    end

    fprintf('RMSE OSE f/cast: %4.4i/%2.2i - %4.4i/%2.2i\n',...
    YRPLT(1,1),YRPLT(1,2),YRPLT(end,1),YRPLT(end,2));
    fprintf('Data extraction: %s - %s\n',...
    datestr(YRPLT(1,3)),datestr(YRPLT(end,3)));

    pthf=sprintf('/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/forecast/%s/%4.4i%2.2i/',...
      esimH,YR,imo);
     
% Loop: forecast 3mo
    iyr=YRPLT(1,1);
    iday=YRPLT(1,2);
    dJ1=datenum(iyr,1,1);
    dnmb0=dJ1+iday-1;

    frmseout = sprintf('%sRMSE%4.4im_OSEfcst_%s_%i%2.2i.mat',pthout,abs(Z0),esimH,iyr,imo);

    clear RMSERR
    RMSERR.Name       = esimH;
    RMSERR.Name_short = esimH;
    RMSERR.Indx_hycom = Iocn;
    RMSERR.HYCOM_dim  = [mh,nh];

    dmean=1;

% Get persistence:
% 1st day f/cast (not saved) = dnmb day from hindcast
    iyr=YRPLT(1,1);
    iday=YRPLT(1,2);
    dJ1=datenum(iyr,1,1);
    dnmb=dJ1+iday-1;
    DV=datevec(dnmb);
    dd=dnmb;
    HR=0;
    pthi1 = sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_%s/',...
                  iyr,esimH);
    fifa_p = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthi1,iyr,iday,HR);
    fifb_p = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthi1,iyr,iday,HR);
    fld='srfhgt';
% Forecast:   
    ssh_prst = sub_getSSH_hycom(fifa_p,fifb_p,INH,dmean);

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
      end
    end

    clear LCXY LCE
    nrc=size(YRPLT,1);
%    tic;
    for ip=1:nrc       % 1 forecast 
      iyr  = YRPLT(ip,1);
      iday = YRPLT(ip,2);
      dJ1  = datenum(iyr,1,1);
      dnmb = dJ1+iday-1;
      DV   = datevec(dnmb);
      mo   = DV(2);
      mday = DV(3);
      dm   = mday;
      dd   = dnmb;
      HR   = 0;

      if f_mat==2
        if ip<=last_rec
          fprintf('OSE fcast %s: %s Record exist skipping %i/%i/%i\n',...
             esimH,datestr(dnmb0),DV(1:3));
          continue
        else
          fprintf('\n\n ==========================================\n')
          fprintf('Continue from ii=%i %s %i/%i/%i\n\n',ip,esimH,DV(1:3));
          f_mat=1;
          last_rec=0;
        end
      end
%      keyboard
      tic;
      fprintf('Reading: %s; esim=%s\n',datestr(dnmb),esimH);

      fifa = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthf,iyr,iday,HR);
      fifb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthf,iyr,iday,HR);

      if ip==1
        fifa = fifa_p;
        fifb = fifb_p;
      end
%      if ~exist(fifa,'file')
%        fprintf('Missing %s\n',fifa);
%        continue
%      end

      fld='srfhgt';
      dmean=1;
% Forecast:   
      ssh_fcst = sub_getSSH_hycom(fifa,fifb,INH,dmean);
%
% Control fields: 50.1 GOMu daily mean SSH
      fprintf('Reading control SSH 50.1 GOMu Reanalysis \n');
      fsshmn = sprintf('%sssh_daily_GOMu501_%4.4i%2.2i.mat',pthmat_sshR,iyr,mo);
      if isempty(SSHM)
        load(fsshmn);
      end
      if SSHM(1).YR ~= iyr | SSHM(1).Mo ~= mo;
        load(fsshmn);
      end
      llc = length(SSHM);
      icR = 0;
      for icR = 1:llc
        if SSHM(icR).day == dm;
          break;
        end
      end
      if SSHM(icR).day ~= dm,
        fprintf('%s not found in daily mean SSH 50.1GOMu\n',datestr(dnmb));
        continue;
      end
      sshR = SSHM(icR).ssh_mn;
      if dmean > 0;
        mnR = nanmean(sshR(INR));
        sshR = sshR - mnR;
      end

% ------------------
% Calc RMSE
% ------------------
      RMSE =[];
      RMSEP = [];
      cnn = 0;
      for inx=1:Nocn
        if mod(inx,10000)==0,
          fprintf(' Calculating RMSE: %4.1f%% done ...\n',inx/Nocn*100);
        end

        ih0=Iocn(inx);
        Iocn_gom = INDX.IAS_HYCOM_Iocn; % ocean indices of GOM for IAS->GOMu mapping
        II=find(Iocn_gom == ih0);
        if isempty(II);
          RMSE(inx)=nan;
          continue;
        end; % pnt outside saved indices

        Igomu = INDX.I_GOMu(II,:);
        Jgomu = INDX.J_GOMu(II,:);

        f_chck=0;
        if f_chck==1
          figure(10); clf;
          hold on;
          plot(LON(ih0),LAT(ih0),'o');

          for itt=1:4
            if itt>1
            plot(LONR(Jgomu(itt),Igomu(itt)),LATR(Jgomu(itt),Igomu(itt)),'k*');
            else
            plot(LONR(Jgomu(itt),Igomu(itt)),LATR(Jgomu(itt),Igomu(itt)),'r*');
            end
          end
        end

  %
  % Take the smallest RMSD
  % out of the surrounding closest GOMu points:
        Err  = [];
        eH   = ssh_fcst(ih0);
        eHp  = ssh_prst(ih0);
        for itt=1:4
          in        = Igomu(itt);
          jn        = Jgomu(itt);
          Err(itt)  = (eH-sshR(jn,in)).^2;
          ErrP(itt) = (eHp-sshR(jn,in)).^2;
        end

        cnn=cnn+1;
        RMSE(inx)  = min(Err);
        RMSEP(inx) = min(ErrP);  

  % Save GOMu ssh on IAS grid:
        imm = find(Err==min(Err),1);
        sshGOMrmp(inx) = sshR(Jgomu(imm),Igomu(imm));

      end;  % ocean points

      rms  = sqrt(nansum(RMSE)/cnn);
      rmsp = sqrt(nansum(RMSEP)/cnn);


% Spatial map of RMSE;
      f_chck=0;
      if ip == f_chck
        fprintf(' Plotting check out ...\n');
        ERR=zeros(m,n)*nan;
        ERR(Iocn)=sqrt(RMSE);
        figure(11); clf;
        pcolor(ERR); shading flat;

        axis('equal');
        set(gca,'xlim',[1 620],...
                'ylim',[390 810]);
        colorbar

        sshGOM = ERR*nan;
        sshGOM(Iocn) = sshGOMrmp;
        figure(12); clf;
        pcolor(sshGOM); shading flat;
        axis('equal');
        set(gca,'xlim',[1 620],...
                'ylim',[390 810]);
        colorbar

  % Difference:
        figure(13); clf;
        pcolor(sqrt((sshGOM-ssh_fcst).^2)); shading flat;
        axis('equal');
        set(gca,'xlim',[1 620],...
                'ylim',[390 810]);
        colorbar

        keyboard
      end

% Note saved squqared error of daily fields -  use after for computing RMSE
% after all days are processed
      RMSE2_fcst = RMSE(:);
      RMSE2_prst = RMSEP(:);
      RMSERR.ERR_squared(:,ip) = RMSE2_fcst;
      RMSERR.ERRprst_squared(:,ip) = RMSE2_prst;  % persistence squared err

%      fprintf('1 record %6.3f min, RMSE=%8.4f \n',toc/60,rms);
      fprintf('1 record %6.3f min, RMSE=%8.4g, PerstRMSE=%8.4g \n',toc/60,rms,rmsp);

      if mod(ip,15)==0 & f_mat==1
        fprintf('Saving %s\n',frmseout);
        save(frmseout,'RMSERR');
      end
    end   %  days in forecast 

    if f_mat==1
      fprintf('Saving %s\n',frmseout);
      save(frmseout,'RMSERR');
    end

  end     %  months with individaul forecasts
end       % years



