% Calculate RMSE between OSE (PIES or noPIES)
% 0.04 HYCOM GOMu reanalysis 50.1
% HYCOM + NCODA Gulf of Mexico 1/25° Reanalysis
% Tidal forcing is included!
% Need to calc daily average SSH: calc_meanSSH_gomu501.m
%
% PIES and  no PIES are real observations)
%
%


addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

ys = 2011;
ye = ys;
esimH = 'PIES';
%esimH = 'noPIES';

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
pthmat_sshR  = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthtopo_rnls = '/nexsan/archive/GOMu0.04_501/topo/';
pthout       = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst     = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';
pthrnls      = '/nexsan/people/abozec/TSIS/data/aviso/GRIDDED/';


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
for iyr = ys:ye
  pthi1 = sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_%s/',...
                  iyr,esimH);
  pthiR = sprintf('/nexsan/archive/GOMu0.04_501/data/netcdf/%4.4i/',iyr);
  fcstname = sprintf('OSE_hindcast');
  frmseout = sprintf('%sRMSE%4.4im_OSEhindcast%s_GOMu_%i.mat',pthout,abs(Z0),esimH,iyr);

  fprintf('\n\n %s\n',frmseout);
  fprintf(' Forecast    : %s %i \n',fcstname,iyr);
  fprintf(' Input data  : %s\n',pthi1);
  fprintf(' Control data: %s\n',pthiR);
  fprintf(' Output saved: %s\n',frmseout);

  clear RMSERR
  RMSERR.Name       = esimH;
  RMSERR.Name_short = esimH;
  RMSERR.Indx_hycom = Iocn;
  RMSERR.HYCOM_dim  = [mh,nh];

  dmean=1;

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
        
  YRPLT=[];
  cc=0;
  id1=1;
  id2=365;
  ic=mod(iyr,4);
  if ic==0 & id2==365, id2=366; end;
  if iyr>ys
    id1=1;
  end
  for iday=id1:id2
    cc=cc+1;
    jd1=datenum(iyr,1,1);
    dnmb=jd1+iday-1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=iday;
    YRPLT(cc,3)=dnmb;
    DV(cc,:) = datevec(dnmb);
  end
  nrc=size(YRPLT,1);

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
    tic;

%keyboard
    yr  = DV(ii,1);
    mo  = DV(ii,2);
    dm  = DV(ii,3);
    dnmb= YRPLT(ii,3);
    iday= YRPLT(ii,2);
    HR  = 0;
    fina = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthi1,iyr,iday,HR);
    finb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthi1,iyr,iday,HR);

    if ~exist(fina,'file')
      fprintf('Missing %s\n',fina);
      continue
    end

    ssh = sub_getSSH_hycom(fina,finb,INH,dmean);

% Control SSH - 50.1 GOMu reanalysis
    fprintf('Reading control SSH 50.1 GOMu Reanalysis \n');
    fsshmn = sprintf('%sssh_daily_GOMu501_%4.4i%2.2i.mat',pthmat_sshR,yr,mo);
    if isempty(SSHM)
      load(fsshmn);
    end
    if SSHM(1).YR ~= yr | SSHM(1).Mo ~= mo;
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
      eH   = ssh(ih0);
      for itt=1:4
        in=Igomu(itt);
        jn=Jgomu(itt);
        Err(itt)=(eH-sshR(jn,in)).^2;
      end

      cnn=cnn+1;
      RMSE(inx)=min(Err);
%      RMSEP(inx)=min(ErrP);  no persistence for hindcast

% Save GOMu ssh on IAS grid:
      imm = find(Err==min(Err),1);
      sshGOMrmp(inx) = sshR(Jgomu(imm),Igomu(imm));

    end;  % ocean points


%    RMSE = (ssh(Iocn)-sshn(Iocn)).^2;  % f/cast - control run
    rms  = sqrt(nansum(RMSE)/cnn);

% Spatial map of RMSE;
    f_chck=0;
    if f_chck==1
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
      pcolor(sqrt((sshGOM-ssh).^2)); shading flat;
      axis('equal');
      set(gca,'xlim',[1 620],...
              'ylim',[390 810]);
      colorbar


      keyboard
    end

    RMSE=RMSE(:);
    RMSERR.ERR_squared(:,ii)=RMSE;
%    RMSERR.ERRprst_squared(:,ii)=RMSEP;  % persistence squared err

    fprintf('1 record %6.3f min, RMSE=%8.4g \n',toc/60,rms);

    if mod(ii,10)==0 & f_mat==1  % f_mat=2 should change to f_mat=1
      fprintf('Saving %s\n',frmseout);
      save(frmseout,'RMSERR');
    end
  end  % year loop

  if f_mat==1
    fprintf('Saving %s\n',frmseout);
    save(frmseout,'RMSERR');
  end

end

