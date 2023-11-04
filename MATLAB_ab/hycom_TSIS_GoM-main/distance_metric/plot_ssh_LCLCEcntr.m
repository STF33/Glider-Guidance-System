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
ixx = 3; % hindcast/free run # expt

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/'; % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
btx='plot_ssh_LCLCEcntr.m';

% HYCOM-TSIS hindcast experiments:
fhnd = 'hycom_tsis_expts.mat';
load(fhnd);
Nruns = length(EXPT);



% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;


YPLT=[];
cc=0;
for iy=2011:2012
 for dd=1:365
  if iy==2011 & dd==1; continue; end;
  if iy==2012 & dd>182,
   break;
  end
  dnmb=datenum(iy,1,1)+dd-1;
  dv=datevec(dnmb);
  cc=cc+1;
  YPLT(cc,1)=iy;
  YPLT(cc,2)=dv(2);
  YPLT(cc,3)=dv(3);
  YPLT(cc,4)=dd;
  YPLT(cc,5)=dnmb;
 end
end

nrc=cc;


% GoM region, NEMO:
GOMN=[     100     365
     651     337
    1091     687
    1246     798
    1512     881
    1787     998
    1904    1292
    1710    1914
     23    1920
      8     748];


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
TMN = LCN(1).TM; % Nemo
nrc = length(LCN(1).XY);

%cmp=flipud(colormap_cold(360));
% Colormap
ncc=100;
cl1=[0,0.3,0.6];
cl2=[1,1,1];
clr1=mix_2colors(cl1,cl2,ncc);

%cl1=cl2;
%cl2=[0.6,1,0.8];
cl1=cl2;
cl2=[0.6,0.9,0.7];
clr2=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[1,0.8,0.2];
clr3=mix_2colors(cl1,cl2,ncc);
cmp = [clr1;clr2;clr3];
cmp = smooth_colormap(cmp,5);



for ii=190:190
  irc=ii;
  tic;

  yr  = YPLT(ii,1);
  mo  = YPLT(ii,2);
  dm  = YPLT(ii,3);
  dyr  = YPLT(ii,4);
  dnmb = YPLT(ii,5);
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
  fprintf('1 fcst 3: %6.2f min\n\n',toc/60);

 %
 % Get contours from NEMO and HYCOM hindcast:
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;

 %  if IEXX(ixx)==0, continue; end;

  fmatout = sprintf('%shycom_LCcontour_%2.2i.mat',pthmat,ixx);
  fprintf('Loading %s\n',fmatout);
  load(fmatout);

  TM = LCXY.TM;

  dm1= TMN(irc);
  ihc= find(TM==dm1);

  xn = LCN(1).XY(irc).X; % NEMO LC contour
  yn = LCN(1).XY(irc).Y;

  [Xc,Yc] = sub_combineLCLCE(xn,yn,LCEN,lnW,irc,10);

  if ~isempty(ihc);
    xh1 = LCXY.XY(ihc).X;
    yh1 = LCXY.XY(ihc).Y;

   [Xhc,Yhc] = sub_combineLCLCE(xh1,yh1,LCE,lnW,irc,10);

    P = [Xc,Yc];
    Q = [Xhc,Yhc];
    mhd1 = modified_hausdorff_distance(P,Q,'geo');
  else
    mhd1 = nan;
    Xhc = [];
    Uhc = [];
  end;

%  MHD(irc,1) = mhd1;


%
% PLOTTING
  figure(1); clf;
  set(gcf,'Position',[1573         560         949         751]);
  axes('Position',[0.08 0.4 0.8 0.5]);
   % hindcast
  pcolor(LONN,LATN,enm); shading flat;
  colormap(cmp);
  clm1=-0.3;
  clm2=0.6;
%  caxis([-0.5 0.5]);
  caxis([clm1, clm2]);
  hold;
  plot(Xc,Yc,'.','Color',[0 0 0]);
  plot(Xhc,Yhc,'.','Color',[1 0 0]);

  contour(LON,LAT,HH,[0 0],'k');

  plot([-96.8 -95.8],[30.2 30.2],'r-','Linewidth',1.5);
  text(-95.4,30.2,sprintf('HYCOM '));
  
  axis('equal');
  set(gca,'xlim',[-97.8 -80.7],...
   'ylim',[18.5 31],...
   'xtick',[-98:2:-82],...
   'ytick',[18:2:32],...
   'tickdir','out',...
   'Fontsize',12);
  hb=colorbar;
  set(hb,'Position',[0.81 0.35 0.015 0.5],...
   'Ticks',[-1.0:0.2:1.0],...
   'Fontsize',12);
  
 % stl=sprintf('ssh, NEMO, LC/LCE contours (17cm) HYCOM hindcast %2.2i, %s',ixx-1,datestr(dnmb));
  stl=sprintf('ssh, NEMO, LC/LCE (17cm) HYCOM hindcast %2.2i, %s, mhd=%5.1f',...
              ixx-1,datestr(dnmb),mhd1);
  title(stl);

  bottom_text(btx,'pwd',1,'Position',[0.08 0.3 0.4 0.04]);
%keyboard
end



  
  
  
  
  



