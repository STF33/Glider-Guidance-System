% Plot RMSE maps for OSE 3-mo f/casts
% with PIES vs no PIES (real observations)
% RMSE computed for OSE f/cast SSH vs 50.1 GOMu reanalysis
% RMSE computed at calc_rmse_ssh_OSEfcst_501GOMu.m
% 

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

%esimH = 'PIES';
esimH = 'noPIES';

%
% Specify depth limit for RMSE calculation
Z0 = -200; % only deep GoM

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';



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


% Combine all f/casts by 1,2,3 months
cnn=0;
for YR=2009:2010
  im1=1;
  im2=12;
  if YR==2009
    im1=5;
  end

  for imo=im1:im2
%    frmseout = sprintf('%sRMSE%4.4im_OSEfcst_%i%2.2i.mat',pthout,abs(Z0),YR,imo);
    frmseout = sprintf('%sRMSE%4.4im_OSEfcst_%s_%i%2.2i.mat',pthout,abs(Z0),esimH,YR,imo);
    fprintf('Loading %s\n',frmseout);
    load(frmseout);

    cnn=cnn+1;
    dmm=sqrt(RMSERR.ERR_squared);

    if cnn==1
      RMSE=dmm;
    else
      [a1,a2]=size(RMSE);
      [d1,d2]=size(dmm);
      if d2<a2  % last day is missing in several f/casts
        for kk=d2+1:a2
          dmm(:,kk)=dmm(:,kk-1);
        end
      end
      RMSE=RMSE+dmm;
    end
  end
end

RMSE=RMSE/cnn;

Hnd_name = RMSERR.Name_short;

btx = 'plot_RMSE_OSEfcst_map.m';


for imo=1:3
% Average by months
  if imo==1
    rmse1=nanmean(RMSE(:,1:31),2);
  elseif imo==2
    rmse1=nanmean(RMSE(:,32:61),2);
  else
    rmse1=nanmean(RMSE(:,62:91),2);
  end

  fprintf('Plotting month %i\n',imo);
  nfg=imo;
  stl = sprintf('OSE %s 3mo f/casts, Z0=%im, Mo=%i',...
               Hnd_name,Z0,imo);

  sub_plot_rmse(nfg,HH,LON,LAT,Iocn,rmse1,stl);
  bottom_text(btx,'pwd',1);
end


