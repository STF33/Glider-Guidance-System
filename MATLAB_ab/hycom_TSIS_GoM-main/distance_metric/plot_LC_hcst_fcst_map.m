% Plot SSH and LC contours for each day
% with MHD scale for selected forecast
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

s_fig=1;

TISL=0.1;   % target isoline for LC identification
esim='PIES';
ys=2009; % 1 year at a time
im1=5;  % forecast start month
%if ys==2009, im1=5; end;



rg=9806;  % convert pressure to depth, m
huge=1e20;


ptht   = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/figLC/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx='plot_LC_hcst_fcst_map.m';

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


mday=[31;28;31;30;31;30;31;31;30;31;30;31];
if mod(ys,4)==0; mday(2)=29; end;

LC = struct; 
nlc=0;  % counter of forecasts 

imo=im1;

YRPLT=[];

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


iyr=YRPLT(1,1);
iday=YRPLT(1,2);
dJ1=datenum(iyr,1,1);
dnmb0=dJ1+iday-1;

nrc=size(YRPLT,1);

cmp=flipud(colormap_cold(360));

cc=0;
for ip=1:nrc
  tic;
  iyr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
  dJ1=datenum(iyr,1,1);
  dnmb=dJ1+iday-1;
  DV=datevec(dnmb);
  dd=dnmb;
%    mo=DV(2); % this is different from imo - 1st month of forecast
  HR=0;

  fprintf('Reading: %s; esim=%s\n',datestr(dnmb),esim);

  pthi=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_%s/',...
	       iyr,esim);
  pthf=sprintf('/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/forecast/%s/%4.4i%2.2i/',...
		 esim,iyr,imo);
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
    LC.MHD(ip)=nan;
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

% MHD
  P=[xx,yy];
  Q=[xf,yf];
  mhd  = modified_hausdorff_distance(P,Q);
  LC.MHD(ip)=mhd;


  clf;
  axes('Position',[0.08 0.4 0.8 0.5]);
    % hindcast
  pcolor(LON,LAT,sshf0); shading flat;
  colormap(cmp);
  caxis([-0.5 0.5]);
  hold;
  plot(xx,yy,'k-','Linewidth',2.2);
  plot(xf,yf,'r-','Linewidth',1.5);

  contour(LON,LAT,HH,[0 0],'k');

  plot([-96.8 -95.8],[30.2 30.2],'r-','Linewidth',1.5);
  text(-95.4,30.2,sprintf('Forecast, %2.2i cm',TISL*100));
  
  axis('equal');
  set(gca,'xlim',[-97.8 -80.7],...
	  'ylim',[18.5 31],...
	  'xtick',[-98:2:-82],...
	  'ytick',[18:2:32],...
	  'tickdir','out',...
	  'Fontsize',12);
  hb=colorbar;
  set(hb,'Position',[0.81 0.43 0.015 0.5],...
	 'Ticks',[-1:0.25:1],...
	 'Fontsize',12);
  
  dfcst=dnmb-dnmb0+1;
  stl=sprintf('ssh, H/cast %s, %s, fcst day %i',esim,datestr(dnmb),dfcst);
  title(stl);

% MHD  
  MHD=LC.MHD;
  MHD(1)=0;
  Td=[1:length(MHD)];
  axes('Position',[0.11 0.13 0.82 0.18]);
  hold on
  plot(Td,MHD,'Linewidth',2);
  stl=sprintf('MHD LC, %s, Fcst vs Hcst, %i/%2.2i',esim,iyr,im1);
  title(stl);
  xlabel('Forecast day');
  yll=ceil(max(MHD)*10)/10+0.02;
  if yll==0, yll=0.5; end;
  set(gca,'tickdir','out',...
	'xlim',[1 93],...
	'ylim',[0 yll],...
	'xtick',[0:5:95],...
	'xgrid','on',...
	'ygrid','on');

%  bottom_text(btx,'pwd',1);
%  keyboard
  
  if s_fig==1
    set(gcf,'Position',[1497 577 1000 754]);
    bottom_text(btx,'pwd',1);
    fgnm=sprintf('%sssh_fcst_hcst%2.2i_%i%2.2i-%3.3i.png',...
		 pthfig,TISL*100,iyr,imo,cc);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r200',fgnm);
  end
  
  
  fprintf('1 fcst 3: %6.2f min\n\n',toc/60);

end % loop for 3-mo forecast






    
    
    
    
    
    
    
    



