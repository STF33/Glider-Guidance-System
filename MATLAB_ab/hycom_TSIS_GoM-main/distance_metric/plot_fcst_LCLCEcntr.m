% Plot SSH and LC LCE contours used for calculating MHD
% derived from HYCOM forecasts
% extr_lc_hycom_fcst.m
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

s_fig=0; % 1 - plot from frame=1, >1 - skip frames <= s_fig

% Select forecast to plot 
iFcst  = 8; % #2, #3, #6, #7, #8 
itime0 = 1; % 1 - May/June 2011, 2 - Jan.Febr 2012
irun0  = 1; % 7 runs each f/cast series

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthfrm = '/Net/kronos/ddmitry/hycom/TSIS/FIG/frames_fcst/'; 

btx='plot_fcst_LCLCEcntr.m';


% Array with NEMO LC
%  fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat); OLD
fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat1);
fprintf('Loading %s\n',fmatout);
BB=load(fmatout);
LCN=BB.LCXY;  % LC contour
LCEN=BB.LCE;  % LCE contours
TMN = LCN(1).TM;  % Nemo

clear LCXY LCE

%
% Info for HYCOM forecasts
if iFcst<10
  FCST = sub_fcst_info(iFcst);
else
  FCST = sub_fcstPrdct_info(iFcst);
end
Nhnd = FCST.Nhind;  % initial cond from the  hindcast #
hnd_name = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
irun1 = FCST.run1;
irun2 = FCST.run2;
ntime=FCST.ntime; % 2 time windows for forecasts
NTime=ntime;

% All hindcast experiments
load('hycom_tsis_expts.mat');

% Array of forecast runs
icc=0;
for it=1:NTime
  for irun=irun1:irun2
    icc=icc+1;
    FRUNS(icc).Nhndcst = Nhnd;
    FRUNS(icc).TimePeriod = it;
    FRUNS(icc).Nfcst = irun;
    FRUNS(icc).matname = sprintf('hycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',Nhnd,it,irun);
    FRUNS(icc).Hind_name = EXPT(Nhnd).Name;
  end
end

%fmat = 'hycom_info_forecasts.mat';
%fprintf('Loading forecast info %s\n',fmat);
%load(fmat);

% Find forecast
Ntot = length(FRUNS);
fcstout = [];
for ixx=1:Ntot
  nh=FRUNS(ixx).Nhndcst;
  it=FRUNS(ixx).TimePeriod;
  nr=FRUNS(ixx).Nfcst;

  if nh==iFcst & it==itime0 & nr==irun0
    flnm = FRUNS(ixx).matname;
    fcstout = sprintf('%s%s',pthmat,flnm);
    break;
  end;
end

if isempty(fcstout)
  fprintf('Could not find forecast %2.2i-%2.2i%2.2i\n',iFcst,itime,irun);
end

fprintf('Loading Forecast contours %s\n',fcstout);
load(fcstout);

TM  = LCXY.TM;
nrc = length(TM);



% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;


% GoM region, NEMO:
GOMN=[         100         365
         651         337
        1091         687
        1246         798
        1512         881
        1787         998
        1904        1292
        1710        1914
          23        1920
           8         748];


fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';

cntr=0;
LONN=[];
LATN=[];
INN=[];
lnW = -190; % cutoff longitude
LCE1 = 10; % when cobine controus, pick only 1 eastern most LCE

cmp=flipud(colormap_cold(360));


figure(1); clf;
%set(gcf,'Position',[1249 478 1158 825]); % opal
set(gcf,'Position',[21 1 1156 801]); % Bookpro

xPs  = [];
yPs  = [];
xhPs = [];
yhPs = [];
%for ii=1:nrc
for ii=20:20
  dm1= TM(ii);
  if ii<=s_fig & s_fig>1; 
    fprintf('%i Skipping day %s\n',ii,datestr(dm1));
    continue;
  end

  tic;
  irc= find(TMN==dm1);
  xn = LCN(1).XY(irc).X; % NEMO LC contour
  yn = LCN(1).XY(irc).Y;

  dv   = datevec(dm1);
  yr   = dv(1);
  mo   = dv(2);
  dm   = dv(3);
  dyr  = dm1-datenum(yr,1,1)+1;
  dnmb = dm1;
  iday = dyr;

  DV = datevec(dnmb);
  dnmb1=datenum(yr,mo,1);
  dnmb2=dnmb1+32;
  v2=datevec(dnmb2);
  dnmb2=datenum(v2(1),v2(2),1);
  d2=dnmb2-datenum(yr,mo,1);

  fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
  fprintf('Reading NEMO: %s\n',fnemo);

  fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);

  if isempty(LONN)
    fmesh=sprintf('%smesh_mask.nc',fpwd);
    dmm = ncread(fmesh,'nav_lon');
    LONN = dmm';
    dmm = squeeze(ncread(fmesh,'nav_lat'));
    LATN = dmm';

    [mm,nn] = size(LONN);

    [XM,YM]=meshgrid([1:nn],[1:mm]);
    INN = inpolygon(XM,YM,GOMN(:,1),GOMN(:,2));
    clear XM YM
  end


  enm = squeeze(ncread(fin,'ssh',[1 1 dm],[Inf Inf 1]));
  enm = enm';
  I=find(enm==0);
  enm(I)=nan;

  % Subtract spatial mean ssh
  dmm=enm;
  dmm(INN==0)=nan;
  sshM=nanmean(nanmean(dmm));
  enm = enm-sshM;


% Save NEMO cntr for day1 of fcst
  if isempty(xPs);
    irc0 = find(TMN==TM(1));
    xn = LCN(1).XY(irc0).X; % NEMO LC contour
    yn = LCN(1).XY(irc0).Y;

    [Xc,Yc] = sub_combineLCLCE(xn,yn,LCEN,lnW,irc0,LCE1);

    xPs=Xc;
    yPs=Yc;
  end

  xn = LCN(1).XY(irc).X; % NEMO LC contour
  yn = LCN(1).XY(irc).Y;

  [Xc,Yc] = sub_combineLCLCE(xn,yn,LCEN,lnW,irc,LCE1);

% Save HYCOM 1st contour for persistence
  ihc=ii;
  if isempty(xhPs)
    xh1 = LCXY.XY(1).X;
    yh1 = LCXY.XY(1).Y;

    [Xhc,Yhc] = sub_combineLCLCE(xh1,yh1,LCE,lnW,1,LCE1);

    xhPs = Xhc;
    yhPs = Yhc;
  end

  xh1 = LCXY.XY(ihc).X;
  yh1 = LCXY.XY(ihc).Y;

  [Xhc,Yhc] = sub_combineLCLCE(xh1,yh1,LCE,lnW,ihc,LCE1);



%
% PLOTTING
  clf;
  axes('Position',[0.08 0.3 0.8 0.6]);
    % hindcast
  pcolor(LONN,LATN,enm); shading flat;
  colormap(cmp);
  caxis([-0.5 0.5]);
  hold;
  plot(Xc,Yc,'.','Color',[0 0 0]);   % NEMO contour
  plot(Xhc,Yhc,'.','Color',[1 0 0.3]); % HYCOM contour
  plot(xhPs,yhPs,'.','Color',[0.9 0.8 0]); % HYCOM perist
%  plot(xPs,yPs,'.','Color',[0.6 0.6 0.6]); % NEMO persist

  contour(LON,LAT,HH,[0 0],'k');

% Legend
  di=0.4;
  ix=-97.5; 
  iy=30.5;  
  plot([ix ix+2*di],[iy iy],'k-','Linewidth',1.5);
  text(ix+3*di,iy,sprintf('NEMO '));
  ix=-97.5; 
  iy=30.2;  
  plot([ix ix+2*di],[iy iy],'r-','Linewidth',1.5);
  text(ix+3*di,iy,sprintf('HYCOM '));
  ix=-94.1; 
  iy=30.5;  
%  plot([ix ix+2*di],[iy iy],'-','Color',[0.6 0.6 0.6],'Linewidth',1.5);
%  text(ix+3*di,iy,sprintf('NEMO prst'));
  ix=-94.1; 
  iy=30.2;  
  plot([ix ix+2*di],[iy iy],'-','Color',[0.9 0.8 0],'Linewidth',1.5);
  text(ix+3*di,iy,sprintf('HYCOM prst'));

  
  axis('equal');
  set(gca,'xlim',[-97.8 -80.7],...
   'ylim',[18.5 31],...
   'xtick',[-98:2:-82],...
   'ytick',[18:2:32],...
   'tickdir','out',...
   'Fontsize',12);

  hb=colorbar;
  set(hb,'Position',[0.81 0.32 0.015 0.46],...
  'Ticks',[-1:0.25:1],...
  'Fontsize',12);
  
  nmfcst = sprintf('%2.2i-%2.2i%2.2i',FRUNS(ixx).Nhndcst,FRUNS(ixx).TimePeriod,FRUNS(ixx).Nfcst);
  stl=sprintf('ssh, NEMO, LC/LCE contours HYCOM fcst %s, %2.2i/%2.2i/%4.4i',nmfcst,dv(3),dv(2),dv(1));
  title(stl);

  bottom_text(btx,'pwd',1,'Position',[0.08 0.22 0.4 0.04]);

%keyboard

  if s_fig>0
    fgnm = sprintf('%sLCcntrs_nemo_hycom_fcst%s_%3.3i.png',pthfrm,nmfcst,ii);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r150',fgnm);
  end

  fprintf('Processed time: %6.2f min\n\n',toc/60);


end



    
    
    
    
    



