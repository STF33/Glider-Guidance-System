% Derive LC/LCE contours in the OSE forecasts (3 months f/casts)
% Initialized from OSE hindcasts: PIES/noPIES
% updated code LC_hcst_fcst.m 
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;
startup;

close all
clear

% ----------------
% Flags
% ---------------
s_fmat=1;

TISL=0.17;  % target isoline for LC identification
%esim='PIES';
esim='noPIES';
%ys=2009; % 1 year at a time
YR1=2009;
YR2=2010;
%im1=5;  % forecast start month
%im2=12;
%if ys==2009, im1=5; end;
Z0 = -200;


fprintf('%s %4.3f m\n',esim,TISL);

rg=9806; % convert pressure to depth, m
huge=1e20;


ptht  = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx='LC_OSEhcst_fcst.m';

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);


% GoM region:
GOM=[366  489
  476  531
  583  560
  576  646
  508  827
  336  848
  204  829
  64  798
  19  746
  16  662
  12  578
  25  455
  71  382
  165  356
  281  400];

[XM,YM]=meshgrid([1:n],[1:m]);
%IN = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
Iocn = find(INH==1 & HH<Z0);
Nocn=length(Iocn);
clear XM YM

IN=INH;
Imean = INH;
Imean(HH>Z0)=0;


LC = struct; 
nlc=0; % counter of forecasts 
for YR=YR1:YR2
  ys=YR;
  mday=[31;28;31;30;31;30;31;31;30;31;30;31];
  if mod(ys,4)==0; mday(2)=29; end;

  im1=1;
  im2=12;
  if ys==2009, im1=5; end;
  
  for imo=im1:im2
%    if (YR==2009 & imo~=6) | (YR>2009); continue; end;  % missed f/cast
%    if (YR==2010 & imo<3) | (YR>2009); continue; end;  % missed f/cast
 % Forecast : 3 months
 % Dates for forecast:
    jd1=datenum(ys,1,1);
    id11=datenum(ys,imo,1)-jd1+1; % start Yr. day
    dnmb1 = datenum(ys,imo,1);
% ????????
%    if ys==2009 & imo==6
%     dnmb1 = datenum(ys,imo,16);
%    end
% ?????
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
      dnmb=dnmb+1;
      cc=cc+1;
      DV=datevec(dnmb);
      jd1=datenum(DV(1),1,1);
      yday=dnmb-jd1+1;
      YRPLT(cc,1)=DV(1);
      YRPLT(cc,2)=yday;
      YRPLT(cc,3)=dnmb;
    end

    fprintf('LC contour: %4.4i/%2.2i - %4.4i/%2.2i\n',...
    YRPLT(1,1),YRPLT(1,2),YRPLT(end,1),YRPLT(end,2));
    fprintf('Data extraction: %s - %s\n',...
    datestr(YRPLT(1,3)),datestr(YRPLT(end,3)));

    pthf=sprintf('/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/forecast/%s/%4.4i%2.2i/',...
      esim,YR,imo);

  % Loop: forecast 3mo
    iyr=YRPLT(1,1);
    iday=YRPLT(1,2);
    dJ1=datenum(iyr,1,1);
    dnmb0=dJ1+iday-1;

    clear LCXY LCE
    nrc=size(YRPLT,1);
    cc=0;
    nlc=nlc+1;
    tic;
    for ip=1:nrc       % 1 forecast 
      iyr=YRPLT(ip,1);
      iday=YRPLT(ip,2);
      dJ1=datenum(iyr,1,1);
      dnmb=dJ1+iday-1;
      DV=datevec(dnmb);
      dd=dnmb;
    %   mo=DV(2); % this is different from imo - 1st month of forecast
      HR=0;

      fprintf('Reading: %s; esim=%s\n',datestr(dnmb),esim);


 % Hindcast LC/LCE extracted in extr_lc_OSEhindcast.m
 %    pthi=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_%s/',...
 %          iyr,esim);
 %    fina = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthi,iyr,iday,HR);
 %    finb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthi,iyr,iday,HR);
      fifa = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthf,iyr,iday,HR);
      fifb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthf,iyr,iday,HR);

%    if ~exist(fina,'file')
%     fprintf('Missing %s\n',fina);
%     continue
%    end
%keyboard
      if ~exist(fifa,'file')
        fprintf('Missing %s\n',fifa);
        continue
      end

      cc=cc+1;

      fld='srfhgt';
      dmean=1;
%    [F,n,m,l] = read_hycom(fina,finb,fld);
% Forecast:   
      ssh = sub_getSSH_hycom(fifa,fifb,Imean,dmean);
%
% Derive LC contour:
      dmm=ssh;
      dmm(INH==0)=nan;
      bsgn = -1;
      f_stop = 0;
      Bisol = TISL;
      LCH1 = identify_LC(LON,LAT,dmm,Bisol,'err_stop',f_stop);

      cntr=cc;
     % HYCOM-TSIS
      LCXY.TM(cntr)   = dnmb;
      LCXY.XY(cntr).X  = LCH1(1).xx;
      LCXY.XY(cntr).Y  = LCH1(1).yy;
     % Save LCEs as well:
      LCE(1).TM(cntr)  = dnmb;
      lcc=length(LCH1);
      LCE(1).NumbLCE(cntr) = lcc;
      LCE(1).XY(cntr).X=[];
      LCE(1).XY(cntr).Y=[];
      if lcc>1
        for ilc=2:lcc
          LCE(ilc-1).XY(cntr).X = LCH1(ilc).xx;
          LCE(ilc-1).XY(cntr).Y = LCH1(ilc).yy;
        end
      end

% Plot
      f_chck = 0;
      if cc==20 & f_chck==1
        figure(10); clf;
        pcolor(LON,LAT,dmm); shading flat;
        hold on;
        contour(LON,LAT,dmm,[Bisol Bisol],'k');

        axis('equal');
        set(gca,'xlim',[-98 -80],...
                'ylim',[17 31]);

        xlc = LCXY.XY(cntr).X;
        ylc = LCXY.XY(cntr).Y;
        dstr = datestr(dnmb);
        plot(xlc,ylc,'r.');

        nlce = length(LCE);
        for jj=1:nlce
          neddies = length(LCE(jj).XY);
          if cntr>neddies; continue; end
          xlce = LCE(jj).XY(cntr).X;
          ylce = LCE(jj).XY(cntr).Y;
          if isempty(xlce), continue; end;
          plot(xlce,ylce,'m.');
        end
        stl = sprintf('%s %s fcst day=%i',esim,dstr,cntr);
        title(stl);
        bottom_text(btx,'pwd',1);

        keyboard
      end

      dfcst=dnmb-dnmb0+1;
      fprintf('Fcst day: %i, cntr=%i, Processed 1 rec, %6.4f min\n',dfcst,cntr);

    end % loop for 1  3-mo forecast

    fprintf('1 fcst 3 month: %6.2f min\n\n',toc/60);

    if s_fmat
      fmat=sprintf('%sLC_distance_OSE%s_fcst_%4.4i%2.2i.mat',...
      pthmat,esim,YR,imo);
      fprintf('saving %s\n',fmat);
      save(fmat,'LCXY','LCE');
    end
  end  % loop for all months
end;  % years


  
  
  
  
  
  
  



