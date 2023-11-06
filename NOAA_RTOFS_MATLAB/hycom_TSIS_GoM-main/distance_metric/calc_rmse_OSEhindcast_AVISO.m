NOT FINISHED 
AVISO is too coarse 
Decided to do validation using 0.04 GOMl reanalysis
% Calculate RMSE between OSE (PIES or noPIES)
% and gridded AVISO 
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

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';
pthaviso = '/nexsan/people/abozec/TSIS/data/aviso/GRIDDED/';


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


YRPLT=[];
cc=0;
for iyr=ys:ye
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
end
nrc=size(YRPLT,1);


use 0.04 Gom reanalysis ssh
/nexsan/archive/GOMu0.04_501/data/netcdf/2009



% Find indices of the GoM region on AVISO grid
dindx = '_20210726';
dnmb  = datenum(2010,10,1);
dv    = datevec(dnmb);
flaviso = sprintf('%sdt_global_allsat_phy_l4_%4.4i%2.2i%2.2i%s.nc',...
                  pthaviso,dv(1:3),dindx);
[i1A,i2A,j1A,j2A] = sub_find_indices_AVISO(flaviso,xt1,xt2,yt1,yt2);


for iyr = ys:ye
  pthi1 = sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_%s/',...
                  iyr,esimH);
  pthi2 = pthaviso;
%  pthi2=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_noPIES/',iyr);

  fcstname = sprintf('OSE_hindcast');
  frmseout = sprintf('%sRMSE%4.4im_OSEhindcastAVISO_%i.mat',pthout,abs(Z0),iyr);

  fprintf('\n\n %s\n',frmseout);
  fprintf(' Forecast: %s %i \n',fcstname,iyr);
  fprintf(' Input data: %s\n',pthi1);
  fprintf(' Input data: %s\n',pthi2);
  fprintf(' Output saved: %s\n',frmseout);

  clear RMSERR
  RMSERR.Name='noPIES';
  RMSERR.Name_short='noPIES';
  RMSERR.Indx_hycom=Iocn;
  RMSERR.HYCOM_dim=[m,n];

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

% Control SSH - AVISO
    fprintf('Reading control SSH AVISO \n');
    ssh_aviso = sub_getSSH_aviso(pthaviso,dnmb,i1A,i2A,j1A,j2A,dmean); 
%    sshn = sub_getSSH_hycom(fina,finb,INH,dmean);

    RMSE =[];
    RMSE = (ssh(Iocn)-sshn(Iocn)).^2;  % f/cast - control run
    rms  = sqrt(nansum(RMSE)/Nocn);

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

    fprintf('1 record %6.3f min, RMSE=%8.4g \n',toc/60,rms);

    if mod(ii,20)==0 & f_mat==1  % f_mat=2 should change to f_mat=1
      fprintf('Saving %s\n',frmseout);
      save(frmseout,'RMSERR');
    end
  end  % year loop

  if f_mat==1
    fprintf('Saving %s\n',frmseout);
    save(frmseout,'RMSERR');
  end

end

