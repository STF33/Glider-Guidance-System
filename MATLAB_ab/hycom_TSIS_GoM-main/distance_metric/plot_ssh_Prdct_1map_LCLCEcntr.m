% NEMO - HYCOM
% Plot only HYCOM ssh and contours from NEMO and HYCOM
%
% Phase 3 - Predictability experiments (with initial state
% from interpolated NEMO+GLORYS)
%  all forecast runs are >=10
%  naming convenction: *_fcstNhnd-TPNR, Nhnd = 10, TP = time
%  period =1 (start in 2011 May - June), NR - run # (01 - starts
%  on day 1 of the time period, then 7-day shift)
%
% Plot SSH and LC LCE contours used for calculating MHD
% see mhd_LCLCEcntr_nemo_hycom.m
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

s_fig=0;

%IEXX = zeros(9,1);  % hindcast free run # expt
%IEXX(2:9) = 1;
ixx = 14;  % hindcast/free run # expt
iFcst = ixx;
% Specify forecast run:
nhnd  = iFcst; %initial conditions, 10 - interpolation
itime = 1;  % =1 2011, =2 2012
irun  = 1;  % f/cast # within the time period - 7-day shifted

% Plot for:
irc1 = 89;
irc2 = 89;

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
pthfrm = '/Net/kronos/ddmitry/hycom/TSIS/FIG/frames_fcst2/';

btx='plot_ssh_Prdct_1map_LCLCEcntr.m';

% HYCOM-TSIS hindcast experiments:
fhnd = 'hycom_tsis_expts.mat';
load(fhnd);
Nruns = length(EXPT);

%
% Get LC/LCE contours from the f/cast
fmatout = sprintf('%shycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',pthmat2,nhnd,itime,irun);
fprintf('Loading %s\n',fmatout);
load(fmatout);
TM  = LCXY.TM;
nrec = length(TM);
LCHY=LCXY;
LCEH=LCE;

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

% 
% HYCOM output f/cast:
if iFcst<10
  FCST = sub_fcst_info(iFcst);
else
  FCST = sub_fcstPrdct_info(iFcst);
end

pthhycom = FCST.TIME0(itime).RUN(irun).pthbin; % output dir

% GoM region HYCOM:
% similar to NEMO GOM domain
% for getting same endpoints of the contours
GOM=[366   489
   476   531
   542   560
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
clear XM YM




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

fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat);

cntr=0;
LONN=[];
LATN=[];
INN=[];
lnW = -90; % cutoff longitude

% Array with NEMO LC
fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat);
fprintf('Loading %s\n',fmatout);
load(fmatout);
LCN=LCXY;  % LC contour
LCEN=LCE;  % LCE contours
TMN = LCN(1).TM;  % Nemo
nrc = length(LCN(1).XY);

cmp=flipud(colormap_cold(360));

figure(1); 
set(gcf,'Position',[1364         502         944         840]);
%for irc=1:nrec
for irc=irc1:irc2
  tic;
  dnmb = TM(irc);
  DV   = datevec(dnmb);
  yr   = DV(1);
  mo   = DV(2);
  dm   = DV(3);
  iday = dnmb-datenum(yr,1,1)+1;


% HYCOM   
% Get HYCOM ssh:
  fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthhycom,yr,iday);
  finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthhycom,yr,iday);
  fprintf('Reading %s\n',fina);

  huge  = 2e20;
  rg    = 9806;

  fld = 'srfhgt';
  [F,nn,mm,ll] = read_hycom(fina,finb,fld);
  F(F>huge)=nan;
  ssh=squeeze(F)./(1e-3*rg);  % ssh m
  dmm=ssh;
  dmm(INH==0)=nan;
%  dmm(HH>-200)=nan;
  sshM=nanmean(nanmean(dmm));
  ssh=ssh-sshM;


%  fprintf('Read %6.2f min\n\n',toc/60);

%
% Get contours from NEMO and HYCOM hindcast:
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;

%  if IEXX(ixx)==0, continue; end;
  inm = find(TMN==dnmb);
  dm1 = TMN(inm);
  ihc = irc;

  xn = LCN.XY(inm).X; % NEMO LC contour
  yn = LCN.XY(inm).Y;

%  LCE1=1;
%  [Xc,Yc] = sub_combineLCLCE(xn,yn,LCEN,lnW,inm,LCE1);

 xh1 = LCHY.XY(ihc).X;
 yh1 = LCHY.XY(ihc).Y;



%  MHD(irc,1)  = mhd1;
  Nlcec = length(LCEN);
  Nlcef = length(LCEH);


%
% PLOTTING
  clf;
  axes('Position',[0.08 0.2 0.8 0.7]);
  pcolor(LON,LAT,ssh); shading flat;
  stl=sprintf('ssh, HYCOM fcst %i  %2.2i/%2.2i/%4.4i',ixx,DV(3),DV(2),DV(1));
  colormap(cmp);
  caxis([-0.5 0.5]);
  hold;

  plot(xn,yn,'.','Color',[0 0 0]);
  for ilce = 1:Nlcec
    xlce = LCEN(ilce).XY(inm).X;
    ylce = LCEN(ilce).XY(inm).Y;
    if isempty(xlce); continue; end
    
    plot(xlce,ylce,'.','Color',[0 0 0]);
  end

  plot(xh1,yh1,'.','Color',[1 0 0]);
  for ilce = 1:Nlcef
    nee = length(LCEH(ilce).XY);
    if irc>nee, continue; end;
    xlce = LCEH(ilce).XY(irc).X;
    ylce = LCEH(ilce).XY(irc).Y;
    if isempty(xlce); continue; end

    plot(xlce,ylce,'.','Color',[1 0 0]);
  end
    

  contour(LON,LAT,HH,[0 0],'k');

  plot([-97.5 -96.5],[30.2 30.2],'r-','Linewidth',2.0);
  text(-96.4,30.2,sprintf('HYCOM '));
  plot([-97.5 -96.5],[29.8 29.8],'k-','Linewidth',2.0);
  text(-96.4,29.6,sprintf('NEMO'));

  
  axis('equal');
  set(gca,'xlim',[-97.8 -80.7],...
   'ylim',[18.5 31],...
   'xtick',[-98:2:-82],...
   'ytick',[18:2:32],...
   'tickdir','out',...
   'Fontsize',12);
  title(stl);

  hb=colorbar('SouthOutside');
  set(hb,'Position',[0.25 0.15 0.5 0.015],...
  'Ticks',[-1:0.25:1],...
  'Fontsize',13);
  
  bottom_text(btx,'pwd',1,'Position',[0.08 0.1 0.4 0.04]);

  fprintf('Processed time: %6.2f min\n\n',toc/60);

end



    
    
    
    
    



