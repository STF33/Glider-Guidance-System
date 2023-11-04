% Plot SSH NEMO vs HYCOM expt
% assimilating NEMO 1/100 free running fields
% and NEMO LC contours
%
% Plot HYCOM and NEMO on 1 figure
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


%expt='ugos0'; %
ix0   = 2; % experiment to plot
s_fig = 0; % =0 - do not save; =1 - overwrite; =2 - start after last saved figure.  



% Hindcasts:
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthfig  = '/Net/kronos/ddmitry/hycom/TSIS/FIG/ssh_hycom_nemo/';


btx = 'plot_ssh_hycom_nemo.m';


ii=0;
EXPT = struct;
% Experiments:
% freerun
ii=ii+1;
EXPT(ii).Name = 'FreeRun';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/freerun/';

% No SSS or SST analysis:
% full field SSH + no pies  
ii=ii+1;
EXPT(ii).Name = '2DSSH noPIES noSSS noSST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_full_nopies_newtsis/';
% AVISO tracks SSH + no pies 
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathsSSH noPIES noSSS noSST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_nopies_newtsis/';
% “one track” SSH + no pies
ii=ii+1;
EXPT(ii).Name = '1swathSSH noPIES noSSS noSST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_onetrack_nopies_newtsis/';
% AVISO tracks SSH + ugos pies (no SSS and no SST)
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathSSH ugosPIES noSSS noSST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_ugos_nosss_nosst_newtsis/';

% No SSS analysis but with SST analysis:
%No SSH track + pies distributed all over the GOM domain (1/30 points)
ii=ii+1;
EXPT(ii).Name = 'noSSH allGoMPIES noSSS SST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_nosla_gompies_nosss_newtsis/';
%AVISO tracks SSH + ugos pies (small area distribution)
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathSSH ugosPIES noSSS SST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_ugos_nosss_newtsis/';
% AVISO tracks SSH + extd pies (bigger area distribution)
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathSSH extendedPIES noSSS SST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_extd_nosss_newtsis/';
% No SSH track + pies distributed all over the GOM domain  (1/60 points)
ii=ii+1;
EXPT(ii).Name = 'NoSSH allGoMPIES noSSS SST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_nosla_gompies60_nosss_newtsis/';

Nruns = ii;

for ii=1:Nruns
  if ii~=ix0;
    fprintf('%i : OFF    %s \n',ii,EXPT(ii).Name);
  else
    fprintf('%i : ON ---> %s \n',ii,EXPT(ii).Name);
  end
end





YPLT=[];
cc=0;
for iy=2011:2012
  for dd=190:190
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


%
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
clear XM YM


pthd = EXPT(ix0).path;

cntr=0;
fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';
for ii=1:1:nrc
  tic;
%
% NEMO data:
%ncdump -h https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc
%fin='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc';
%ncid = netcdf.open(fin);
%vid  = netcdf.inqVarID(ncid,'ssh');
%ssh  = netcdf.getVar(ncid,vid,'single');
%[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,vid)
%varname = netcdf.inqVar(ncid,vid)

  yr   = YPLT(ii,1);
  mo   = YPLT(ii,2);
  dm   = YPLT(ii,3);
  dyr  = YPLT(ii,4);
  dnmb = YPLT(ii,5);
  iday = dyr;

  dnmb1=datenum(yr,mo,1);
  dnmb2=dnmb1+32;
  v2=datevec(dnmb2);
  dnmb2=datenum(v2(1),v2(2),1);
  d2=dnmb2-datenum(yr,mo,1);

  cntr=cntr+1;
  foutH=sprintf('%snemo_hycomtsis%s_ssh_%4.4i.png',pthfig,ix0,cntr);


% Skip if exist saved figures 
  if s_fig==2 & exist(foutH,'file'), 
    fprintf('%s %s exist, skipping ...\n',foutN,foutH);
    continue; 
  end;

  fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
  fprintf('Reading NEMO: %s\n',fnemo);

  fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);

  if ~exist('LONN','var')
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

  Bisol=0.17;
  dmm=enm;
  dmm(INN==0)=nan;
%  LCN = identify_LC(LONN,LATN,dmm,Bisol);

% Read in HYCOM ssh
  sday=sprintf('%3.3i',iday);
  hr=0;
  fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd,yr,iday);
  finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd,yr,iday);
  fin=fina;

  ie = exist(fin,'file');

  if ~ie
    fprintf('Missing: %s\n',fin);
    continue;
  end

  fprintf('Reading %s\n',fina);
  fld = 'srfhgt';
  [F,nn,mm,ll] = read_hycom(fina,finb,fld);
  F(F>huge)=nan;
  ssh=squeeze(F)./(1e-3*rg);  % ssh m
%
% Subtract anomaly:
  dmm=ssh;
  dmm(INH==0)=nan;
%  dmm(HH>-200)=nan;
  sshM=nanmean(nanmean(dmm));
  ssh=ssh-sshM;

%
% Derive LC contour:
% 
  dmm=ssh;
  dmm(INH==0)=nan; 
%  LCH = identify_LC(LON,LAT,dmm,Bisol);

  fprintf('Plotting ...\n');
% =====================
% Plotting
% =====================
%  cntr=cntr+1;
% NEMO:
  date_str=sprintf('%4.4i/%2.2i/%2.2i ',yr,mo,dm);

  stt=sprintf('NEMO SSH %s',date_str);
%  if cntr==1
%    figure('Position',[1064 229 1421 1113]); clf;
%  end
    clr1=[0 0.6 1];
    clr2=[0 1 0.8];

    figure(1); clf;
 
    POS = [0.05 0.08 0.44 0.85; ...
           0.53 0.08  0.44 0.85];
    fnb=1;
    c1=-0.5;
    c2=0.5;
    hhN=[];
    xl1=-98;
    xl2=-78.4;
    yl1=14.;
    yl2=31.;
    clb=1;
    ps1=POS(1,:);
    sub_plot_ssh3(ps1,enm,LONN,LATN,hhN,c1,c2,stt,INN,xl1,xl2,yl1,yl2,clb);
%    plot(LCN(1).xx,LCN(1).yy,'-','Color',clr1,'Linewidth',1.5);
%    plot(LCH(1).xx,LCH(1).yy,'-','Color',clr2,'Linewidth',1.5);

  % Legend
%    axes('Position',[0.08 0.22 0.15 0.1]);
%    hold on;
%    plot([0 0.2],[0.1 0.1],'-','Color',clr1,'Linewidth',1.5);
%    text(0.22, 0.1,'NEMO LC','Fontsize',12);
%    plot([0 0.2],[0.3 0.3],'-','Color',clr2,'Linewidth',1.5);
%    text(0.22, 0.3,'HYCOM LC','Fontsize',12);
%    set(gca,'xlim',[-0.05 0.8],'ylim',[0 0.4],...
%            'xtick',[],'ytick',[]);


% ==============
% HYCOM-TSIS
% ==============
    ps1=POS(2,:);
    stt2=sprintf('HYCOM-TSIS SSH %s',date_str);
    clb=0;    

    sub_plot_ssh3(ps1,ssh,LON,LAT,hhN,c1,c2,stt2,INH,xl1,xl2,yl1,yl2,clb);
%    plot(LCN(1).xx,LCN(1).yy,'-','Color',clr1,'Linewidth',1.5);
%    plot(LCH(1).xx,LCH(1).yy,'-','Color',clr2,'Linewidth',1.5);
    nmexpt = EXPT(ix0).Name;
    text(-96.5,15,nmexpt,'Fontsize',14,'Color',[1 0.9 0]);


    set(gcf,'Position',[1000 481 1383 835]);

    bottom_text(btx,'pwd',1,'Position',[0.08 0.05 0.4 0.05]);

    if s_fig>=1
      fprintf('Saving %s\n',foutH);
      print('-dpng','-r150',foutH);
    end
keyboard

  fprintf('Processed 1 rec, %6.4f min\n\n',toc/60);

%keyboard

end

