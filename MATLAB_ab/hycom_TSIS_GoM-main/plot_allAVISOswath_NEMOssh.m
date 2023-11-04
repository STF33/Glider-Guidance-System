% Plot altimetry swaths (all satellites) for 1 day of analysis
% for specific dates
% in the GoM
% show NEMO ssh in the swath 
% For 2011/2012 
%
% In TSIS, we used +/- 4 days of all satellites in the GoM to 
% increase the coverage
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

%pthdat = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/obs/aviso/'; % for 2009
%pthdat  = '/nexsan/people/abozec/TSIS/data/aviso/';  % <-- all data here, gzip
pthdat = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/data/aviso/';  % copied from Alex's dir
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';

% In TSIS - combine 4 days of altimetry to get more data
%dnmb = datenum(2011,06,01);
dnmb0  = datenum(2011,6,1);  % analysis date:
dnmb1  = dnmb0-4;
dnmb2  = dnmb0+4;

% For 2011-2012
% c2 - Cryosat
%  enn - Envisat
% j1g - Jason 1 Geodetic Phase Global Ocean
% j1n - Jason 1 interleaved
% j2 - Jason 2
SAT{1}='c2';
SAT{2}='enn';
SAT{3}='j1g';
SAT{4}='j1n';
SAT{5}='j2';
nst = length(SAT);

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

ln1=min(min(LON));
ln2=max(max(LON));
lt1=min(min(LAT));
lt2=max(max(LAT));

TRCK = struct;
TRCK(1).ST(1).Name = 'Cryosat';
TRCK(1).ST(2).Name = 'Envisat';
TRCK(1).ST(3).Name = 'Jsn1Geod';
TRCK(1).ST(4).Name = 'Jsn1Intrlv';
TRCK(1).ST(5).Name = 'Jason2';

TM = [dnmb1:dnmb2];
nt = length(TM);

for id=1:nt
  dnmb=TM(id);
  dv=datevec(dnmb);

  TRCK(id).TimeMat = dnmb;
  TRCK(id).TimeStr = datestr(dnmb);

 for ist=1:nst
  sat=SAT{ist};

  fnm=sprintf('%s%s/dt_global_%s_phy_vfec_l3_%i%2.2i%2.2i.nc',pthdat,sat,sat,dv(1:3));
  fprintf('Reading %s\n',fnm);
  if ~exist(fnm,'file')
   fprintf('%s does not exist ...\n',fnm);
   continue;
  end

  lat = nc_varget(fnm,'latitude');
  lon = nc_varget(fnm,'longitude');
  J=find(lon>180);
  lon(J)=lon(J)-360;
  time = nc_varget(fnm,'time');

  I=find(lat>lt1 & lat<lt2 & lon>ln1 & lon<ln2);

  TRCK(id).ST(ist).X = lon(I);
  TRCK(id).ST(ist).Y = lat(I);

 end
end

% Plot tracks:
CLR = [0 0.5 0.9; ...
   0.9 0.4 0; ...
   0   1 0.7; ...
   0.8 0 1; ...
   0.4 0.6 0.2];

nl = length(TRCK);

btx='plot_allAVISOswath_NEMOssh.m';

f_chck = 0;
if f_chck==1
%  dnmb = datenum(2011,06,01); % plot 1 day
  figure(1); clf;
  set(gcf,'Position',[1491         726         791         583]);
  axes('Position',[0.08 0.2 0.85 0.72]);
  hold on;
  contour(LON,LAT,HH,[0 0],'color','k','linewidth',2);

  d1 = dnmb1;
  d2 = dnmb2;
  for dnmb=d1:d2
    for id=1:nl
      if TRCK(id).TimeMat == dnmb; break; end;
    end
    dv=datevec(dnmb);
    fprintf('Plotting %i/%i/%i\n',dv(1:3));

    for ist=1:nst
      if ist>length(TRCK(id).ST); break; end; % missing satellites
      x=TRCK(id).ST(ist).X;
      y=TRCK(id).ST(ist).Y;
      clr = CLR(ist,:);
      plot(x,y,'.','Markersize',12,'Color',clr);
    end

  end

  axis('equal');
  set(gca,'tickdir','out',...
      'xlim',[ln1 -78],...
      'ylim',[17 31.9]);
% set(gca,'tickdir','out',...
%     'xlim',[ln1 ln2],...
%     'ylim',[lt1 31.9]);
  dv1 = datevec(d1);
  dv2 = datevec(d2);
  stl=sprintf('%2.2i/%2.2i/%4.4i-%2.2i/%2.2i/%4.4i',dv1(3:-1:1),dv2(3:-1:1));
 %  title(stl,'Fontsize',12);
  text(-97,30.6,stl,'Fontsize',12);

  title('Satellite Swaths, AVISO L3');

  axes('Position',[0.08 0.08 0.25 0.18]);
  hold on;
  for ist=1:nst
    snm = TRCK(1).ST(ist).Name;
    plot([0.02 0.15],[ist/10 ist/10],'-','Linewidth',2,'Color',CLR(ist,:));
    text(0.3,ist/10,snm);
  end
  set(gca,'xtick',[],...
      'ytick',[],...
      'xlim',[0 0.55],...
      'ylim',[0 0.55]);


  bottom_text(btx,'pwd',1);

end

% Read NEMO ssh:
sshN = sub_get_sshNEMO(dnmb0);

% Get ssh values into swath coordinates
lon1 = ln1;
lon2 = -78;
lat1 = 17;
lat2 = lt2;

for id=1:nl
  fprintf('Extracting ssh along the tracks, day id=%i\n',id);
  for ist=1:nst
    fprintf(' ... %s\n',TRCK(1).ST(ist).Name); 
    if ist>length(TRCK(id).ST); break; end; % missing satellites
    x=TRCK(id).ST(ist).X;
    y=TRCK(id).ST(ist).Y;
    if isempty(x); continue; end;

    AA = sub_ssh2swath(x,y,sshN,LON,LAT,lon1,lon2,lat1,lat2);

    TRCK(id).ST(ist).lonH = AA.lonH;
    TRCK(id).ST(ist).latH = AA.latH;
    TRCK(id).ST(ist).Iindx = AA.Iindx;
    TRCK(id).ST(ist).Jindx = AA.Jindx;
    TRCK(id).ST(ist).ssh   = AA.ssh;
  end
end


LMSK = HH*0;
LMSK(HH<0)=1;
cmpH=[0,0,0; 0.8, 0.8, 0.8];

% Colormap for SSH field:
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
cmp_ssh = [clr1;clr2;clr3];
cmp_ssh = smooth_colormap(cmp_ssh,5);


figure(2); clf;
set(gcf,'Position',[1491         726         791         583]);
axes('Position',[0.08 0.2 0.85 0.72]);
hold on;

pcolor(LON,LAT,LMSK); shading flat;
colormap(cmpH);

d1 = dnmb1;
d2 = dnmb2;
for dnmb=d1:d2
  for id=1:nl
    if TRCK(id).TimeMat == dnmb; break; end;
  end
  dv=datevec(dnmb);
  fprintf('Plotting %i/%i/%i\n',dv(1:3));

  for ist=1:nst
    if ist>length(TRCK(id).ST); break; end; % missing satellites
    x   = TRCK(id).ST(ist).lonH;
    y   = TRCK(id).ST(ist).latH;
    ssh = TRCK(id).ST(ist).ssh;
    if isempty(x); continue; end;

% Get rid of land
    I = find(~isnan(ssh));
    x = x(I);
    y = y(I);
    ssh = ssh(I);
    clm1=-0.3;
    clm2=0.6;
    CC  = sub_match_color(ssh,cmp_ssh,clm1,clm2);
    scatter(x,y,18,CC,'filled');
  end
end

axis('equal');
set(gca,'tickdir','out',...
    'xlim',[ln1 -80],...
    'ylim',[18 31.2]);
sttl = sprintf('Sat.Swaths, all AVISO L3, +/-4 days, analysis date=%s',...
     datestr(dnmb0));

% Add UGOS PIES (=1) or UGOS extended PIES (=2)
pthosse = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/OSSE/';
YR=2011;
f_pies = 2;
if f_pies==1
  finp=sprintf('%sugos_mooring_%i.mat',pthosse,YR);
  fprintf('Loading %s\n',finp);
  A=load(finp);
  plat=A.ugos_lat;
  plon=A.ugos_lon;
  sttl = sprintf('All AVISO, PIES, analysis date=%s',...
     datestr(dnmb0));
elseif f_pies==2
  finp=sprintf('%sextd_mooring_2011.mat',pthosse);
  A=load(finp);
  plat=A.extd_lat;
  plon=A.extd_lon;
  sttl = sprintf('All AVISO, extd PIES, analysis date=%s',...
     datestr(dnmb0));
end
if f_pies>0
  plot(plon,plat,'^','Color',[0 0 0],'linewidth',1.5);
end

dv1 = datevec(d1);
dv2 = datevec(d2);
title(sttl);
bottom_text(btx,'pwd',1);


%
% Plot colorbar
% and SSH
%
figure(3); clf;
set(gcf,'Position',[1491         726         791         583]);
axes('Position',[0.08 0.2 0.85 0.72]);
hold on;

pcolor(LON,LAT,sshN); shading flat;
caxis([clm1,clm2]);

contour(LON,LAT,sshN,[0:0.1:0.6],'k-','linewidth',1.2);

axis('equal');
set(gca,'tickdir','out',...
    'xlim',[ln1 -80],...
    'ylim',[18 31.2]);
set(gca,'Color',[0 0 0]);
sttl = sprintf('NEMO SSH, %s\n',datestr(dnmb0));

clb=colorbar('SouthOutside');
colormap(cmp_ssh);
set(clb,'Position',[0.15 0.1 0.66 0.025],...
        'Fontsize',13,...
        'Ticks',[-1:0.1:1],...
        'Ticklength',0.025);

title(sttl);
bottom_text(btx,'pwd',1);


%
%   Plot 1 satellite 
%
figure(4); clf;
set(gcf,'Position',[1491         726         791         583]);
axes('Position',[0.08 0.2 0.85 0.72]);
hold on;

pcolor(LON,LAT,LMSK); shading flat;
colormap(cmpH);

ist0 = 2;
d1 = dnmb1;
d2 = dnmb2;
for dnmb=d1:d2
  for id=1:nl
    if TRCK(id).TimeMat == dnmb; break; end;
  end
  dv=datevec(dnmb);
  fprintf('Plotting %i/%i/%i\n',dv(1:3));

  for ist=ist0:ist0
    if ist>length(TRCK(id).ST); break; end; % missing satellites
    x   = TRCK(id).ST(ist).lonH;
    y   = TRCK(id).ST(ist).latH;
    ssh = TRCK(id).ST(ist).ssh;
    if isempty(x); continue; end;

% Get rid of land
    I = find(~isnan(ssh));
    x = x(I);
    y = y(I);
    ssh = ssh(I);
    clm1=-0.3;
    clm2=0.6;
    CC  = sub_match_color(ssh,cmp_ssh,clm1,clm2);
    scatter(x,y,18,CC,'filled');
  end
end

axis('equal');
set(gca,'tickdir','out',...
    'xlim',[ln1 -80],...
    'ylim',[18 31.2]);
dv1 = datevec(d1);
dv2 = datevec(d2);
title(sprintf('Sat.Swaths, %s AVISO L3, +/-4 days, analysis date=%s',...
      TRCK(1).ST(ist0).Name,datestr(dnmb0)));
bottom_text(btx,'pwd',1);



% 
% Plot Tracks + SST field + UGOS PIES
%
% Colormap for SST field:
% Colormap
ncc=50;
cl1=[0.8,0.2,0.8];
cl2=[0.6,0.2,0.9];
clr1=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[0.,0.4,0.9];
clr2=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[0,0.9,1];
clr3=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[0,0.7,0.4];
clr4=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[1,1,0.3];
clr5=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[0.9,0.4,0];
clr6=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[.9,0.1,0];
clr7=mix_2colors(cl1,cl2,ncc);

cmp_sst = [clr1;clr2;clr3;clr4;clr5;clr6;clr7];
cmp_sst = smooth_colormap(cmp_sst,5);



pthnemo = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
fprintf('Loading NEMO grid %s\n',fgrd);
load(fgrd);


% Read SST:
fldnm = 'toce';
iz0=1;
TT = sub_get_NEMO_TS(dnmb0,fldnm,iz0);

figure(5); clf;
set(gcf,'Position',[1491         726         791         583]);
axes('Position',[0.08 0.2 0.85 0.72]);
hold on;

pcolor(LONN,LATN,TT); shading flat;
colormap(cmp_sst);
ct1=25;
ct2=29;
caxis([ct1 ct2]);


d1 = dnmb1;
d2 = dnmb2;
for dnmb=d1:d2
  for id=1:nl
    if TRCK(id).TimeMat == dnmb; break; end;
  end
  dv=datevec(dnmb);
  fprintf('Plotting %i/%i/%i\n',dv(1:3));

  for ist=1:nst
    if ist>length(TRCK(id).ST); break; end; % missing satellites
    x   = TRCK(id).ST(ist).lonH;
    y   = TRCK(id).ST(ist).latH;
    ssh = TRCK(id).ST(ist).ssh;
    if isempty(x); continue; end;

% Get rid of land
    I = find(~isnan(ssh));
    x = x(I);
    y = y(I);
    ssh = ssh(I);
    clm1=-0.3;
    clm2=0.6;
    CC  = sub_match_color(ssh,cmp_ssh,clm1,clm2);
    scatter(x,y,18,CC,'filled');
  end
end

axis('equal');
set(gca,'tickdir','out',...
    'xlim',[ln1 -80],...
    'ylim',[18 31.2]);
set(gca,'Color',[0 0 0]);
sttl = sprintf('SST, all AVISO, +/-4 days, analysis date=%s',...
     datestr(dnmb0));

% Add UGOS PIES (=1) or UGOS extended PIES (=2)
pthosse = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/OSSE/';
YR=2011;
f_pies = 2;
if f_pies==1
  finp=sprintf('%sugos_mooring_%i.mat',pthosse,YR);
  fprintf('Loading %s\n',finp);
  A=load(finp);
  plat=A.ugos_lat;
  plon=A.ugos_lon;
  sttl = sprintf('SST, All AVISO, PIES, analysis date=%s',...
     datestr(dnmb0));
elseif f_pies==2
  finp=sprintf('%sextd_mooring_2011.mat',pthosse);
  A=load(finp);
  plat=A.extd_lat;
  plon=A.extd_lon;
  sttl = sprintf('SST, All AVISO, extd PIES, analysis date=%s',...
     datestr(dnmb0));
end
if f_pies>0
  plot(plon,plat,'^','Color',[0 0 0],'linewidth',1.5);
end

dv1 = datevec(d1);
dv2 = datevec(d2);
title(sttl);
bottom_text(btx,'pwd',1);

%
% Plot colorbar
%
figure(6); clf;
axes('Position',[0.08 0.5 0.85 0.4]);
pcolor(TT); shading flat
clb=colorbar('SouthOutside');
colormap(cmp_sst);
caxis([ct1 ct2]);

set(clb,'Position',[0.1 0.3 0.66 0.025],...
        'Fontsize',13,...
        'Ticks',[10:1:40],...
        'Ticklength',0.025);



