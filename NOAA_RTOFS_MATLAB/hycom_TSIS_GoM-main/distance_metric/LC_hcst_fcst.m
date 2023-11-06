% Plot distance metric between LC contour
% in hindcast vs forecast
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

TISL=0.10;   % target isoline for LC identification
%esim='PIES';
esim='noPIES';
%ys=2009; % 1 year at a time
YR1=2009;
YR2=2010;
%im1=5;  % forecast start month
%im2=12;
%if ys==2009, im1=5; end;


fprintf('%s %4.3f m\n',esim,TISL);

rg=9806;  % convert pressure to depth, m
huge=1e20;


ptht   = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx='LC_hcst_fcst.m';

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);


% GoM region:
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
IN = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
clear XM YM



LC = struct; 
nlc=0;  % counter of forecasts 
for YR=YR1:YR2
  ys=YR;
  mday=[31;28;31;30;31;30;31;31;30;31;30;31];
  if mod(ys,4)==0; mday(2)=29; end;

  im1=1;
  im2=12;
  if ys==2009, im1=5; end;
  
  for imo=im1:im2
  % Forecast : 3 months
  % Dates for forecast:
    jd1=datenum(ys,1,1);
    id11=datenum(ys,imo,1)-jd1+1; % start Yr. day
    dnmb1 = datenum(ys,imo,1);
    if ys==2009 & imo==6
      dnmb1 = datenum(ys,imo,16);
    end
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

    nrc=size(YRPLT,1);
    cc=0;
    nlc=nlc+1;
    for ip=1:nrc
      iyr=YRPLT(ip,1);
      iday=YRPLT(ip,2);
      dJ1=datenum(iyr,1,1);
      dnmb=dJ1+iday-1;
      DV=datevec(dnmb);
      dd=dnmb;
  %    mo=DV(2); % this is different from imo - 1st month of forecast
      HR=0;

      tic;
      fprintf('Reading: %s; esim=%s\n',datestr(dnmb),esim);

      pthi=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_%s/',...
		   iyr,esim);
      fina = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthi,iyr,iday,HR);
      finb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthi,iyr,iday,HR);
      fifa = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthf,iyr,iday,HR);
      fifb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthf,iyr,iday,HR);

      if ~exist(fina,'file')
				fprintf('Missing %s\n',fina);
				continue
      end
      if ~exist(fifa,'file')
				fprintf('Missing %s\n',fifa);
				continue
      end

      cc=cc+1;

      fld='srfhgt';
      [F,n,m,l] = read_hycom(fina,finb,fld);
      F(F>huge)=nan;
      SSH=squeeze(F)./(1e-3*rg);  % ssh m
  % Delete fields outside GoM
      SSH(IN==0)=nan;
      hMi=nanmean(nanmean(SSH));
      SSH=SSH-hMi;
      sshf0=SSH;
      SSH(HH>-20)=nan;
      LCLCE = identify_LC(LON,LAT,SSH,TISL);
      xx=LCLCE(1).xx;
      yy=LCLCE(1).yy;
      if isempty(xx), error('Could not locate the LC'); end;
      LC(nlc).hcst(cc).xx=xx;
      LC(nlc).hcst(cc).yy=yy;
      LC(nlc).TM(cc)=dnmb;

  % Forecast:    
      [F,n,m,l] = read_hycom(fifa,fifb,fld);
      F(F>huge)=nan;
      SSHf=squeeze(F)./(1e-3*rg);  % ssh m
      SSHf(IN==0)=nan;
      hMf=nanmean(nanmean(SSHf));
      SSHf=SSHf-hMf;
      sshi0=SSHf;
      SSHf(HH>-20)=nan;
      LCLCE = identify_LC(LON,LAT,SSHf,TISL);
      xf=LCLCE(1).xx;
      yf=LCLCE(1).yy;
      if isempty(xx), error('Could not locate the LC f/cast'); end;
      LC(nlc).fcst(cc).xx=xf;
      LC(nlc).fcst(cc).yy=yf;

  % Plot
      f_plt=1;
      if f_plt==1
				% hindcast
				xx=LC(nlc).hcst(cc).xx;
				yy=LC(nlc).hcst(cc).yy;

				% forecast
				xf=LC(nlc).fcst(cc).xx;
				yf=LC(nlc).fcst(cc).yy;

				figure(10); clf;
				axes('Position',[0.08 0.3 0.8 0.6]);
				pcolor(LON,LAT,sshf0); shading flat;
				hold;
				plot(xx,yy,'k-');

				plot(xf,yf,'r--');

% Plot other LCEs
        nlce=length(LCLCE);
        for ilce=2:nlce
          xlce = LCLCE(ilce).xx;
          ylce = LCLCE(ilce).yy;
          plot(xlce,ylce,'b-');
        end

				contour(LON,LAT,HH,[0 0],'k');
				caxis([-0.3 0.3]);
				axis('equal');
				set(gca,'xlim',[-97.9 -80.4],...
					'ylim',[18 31]);
				colorbar

				dfcst=dnmb-dnmb0+1;
				stl=sprintf('ssh, H/cast %s, %s, fcst day %i',esim,datestr(dnmb),dfcst);
				title(stl);

        bottom_text(btx,'pwd',1,'position',[0.1 0.21 0.6 0.05]);

        keyboard

      end

  % 
  % Calculate distance metric
      P=[xx,yy];
      Q=[xf,yf];
      mhd  = modified_hausdorff_distance(P,Q);
      LC(nlc).MHD(cc)=mhd;
      
% Calculate persistence MHD
      if cc==1
				xp0=xx;
				yp0=yy;
      end
      P=[xx,yy];
      Q=[xp0,yp0];
      mhd0  = modified_hausdorff_distance(P,Q);
      LC(nlc).MHD_prst(cc)=mhd0;
      

      dfcst=dnmb-dnmb0+1;
      fprintf('Fcst day: %i, cc=%i, mhd=%8.4f persist mhd=%8.4f\n\n',...
	      dfcst,cc,mhd,mhd0);

    end % loop for 3-mo forecast

    fprintf('1 fcst 3 month: %6.2f min\n\n',toc/60);

    if s_fmat
      fmat=sprintf('%sLC_distance_hcst_fcst_%s_%3.3icm_%i-%i.mat',...
		   pthmat,esim,round(TISL*100),YR1,YR2);
      fprintf('saving %s\n',fmat);
      save(fmat,'LC');
    end
  
  end   % loop for all months
end;   % years


    
figure(1); clf;
axes('Position',[0.08 0.6 0.8 0.3]);
hold on

nL=length(LC);

for ik=1:nL
  MHD=LC(ik).MHD;
  Td=[1:length(MHD)];
  plot(Td,MHD);
end

stl=sprintf('MHD LC, %s, Fcst vs Hcst, %i, %i-%i',esim,iyr,im1,im2);
title(stl);
xlabel('Forecast day');
set(gca,'tickdir','out',...
	'xgrid','on',...
	'ygrid','on');

bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.4 0.05]);
    
    
    
    
    
    
    
    



