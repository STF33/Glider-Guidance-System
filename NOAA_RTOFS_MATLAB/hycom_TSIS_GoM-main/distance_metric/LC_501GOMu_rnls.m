% Derive LC/LCE contours from the 501 GOMu reanalysis 
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
%ys=2009; % 1 year at a time
YR1=2009;
YR2=2011;
Z0 = -200;


ptht    = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat  = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthmat_sshR  = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthfig  = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx     = 'LC_501GOMu_rnls.m';

fmat=sprintf('%sLCLCEcntr_501GOMu_rnls.mat',pthmat);


%
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


SSHM  = [];
YRPLT = [];
cc=0;
for iyr=YR1:YR2
  id1=1;
  id2=365;
  ic=mod(iyr,4);
  if ic==0 & id2==365, id2=366; end;
  if iyr==2009
    id1=datenum(2009,5,1)-datenum(2009,1,1)+1;
  end
  if iyr==2011
    id2 = datenum(2011,4,1)-datenum(2011,1,1)+1;
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



LC = struct; 
nlc=0; % counter of forecasts 
for ii=1:nrc
  tic;

  yr  = DV(ii,1);
  mo  = DV(ii,2);
  dm  = DV(ii,3);
  dnmb= YRPLT(ii,3);
  iday= YRPLT(ii,2);

  pthiR = sprintf('/nexsan/archive/GOMu0.04_501/data/netcdf/%4.4i/',yr);
  fprintf('Reading SSH 50.1 GOMu Reanalysis %i/%2.2i/%2.2i\n',yr,mo,dm);
  fsshmn = sprintf('%sssh_daily_GOMu501_%4.4i%2.2i.mat',pthmat_sshR,yr,mo);
  if isempty(SSHM)
    fprintf('Loading %s\n',fsshmn);
    load(fsshmn);
  end
  if SSHM(1).YR ~= yr | SSHM(1).Mo ~= mo;
    fprintf('Loading %s\n',fsshmn);
    load(fsshmn);
  end
  llc = length(SSHM);
  icR = 0;
  for icR = 1:llc
    if SSHM(icR).day == dm;
      break;
    end
  end
  ssh = SSHM(icR).ssh_mn;
% Subtract anomaly:
  dmm = ssh;
  dmm(INR==0) = nan;
  sshM = nanmean(nanmean(dmm));
  ssh  = ssh - sshM;

% Derive LC contour:
  dmm=ssh;
  dmm(INR==0)=nan;
  bsgn = -1;
  f_stop = 0;
  Bisol = TISL;
  LCH1 = identify_LC(LONR,LATR,dmm,Bisol,'err_stop',f_stop);

  cntr=ii;
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
  if ii==2 & f_chck==1
    figure(10); clf;
    pcolor(LONR,LATR,dmm); shading flat;
    hold on;
    contour(LONR,LATR,dmm,[Bisol Bisol],'k');

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
    stl = sprintf('50.1GOMu Reanalysis %s fcst day=%i',dstr,cntr);
    title(stl);
    bottom_text(btx,'pwd',1);

    keyboard
  end


  ichck=10;
  if mod(ii,ichck)==0
    fprintf('%i rec : %6.2f min\n\n',ichck,toc/60);

    if s_fmat
      fprintf('saving %s\n',fmat);
      save(fmat,'LCXY','LCE');
    end
  end
end;  % years

if s_fmat
  fprintf('saving %s\n',fmat);
  save(fmat,'LCXY','LCE');
end


  
  
  
  
  
  
  



