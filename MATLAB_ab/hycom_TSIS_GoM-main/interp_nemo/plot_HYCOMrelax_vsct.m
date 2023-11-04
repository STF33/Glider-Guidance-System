% Plot vertical sections
% in terpolated fields
% from HYCOM relax files (relax.csh interpolates
% T/S fields from Z-level "climatology" files
% for creating "cold start" state) 
% created from NEMO-GLORYS interpolated
% into IAS HYCOM-TSIS grid and NEMO vertical Z-levels
% 
% Relax files - interpolated into IAS HYCOM hybrid layers
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


fid = [];

fldnm = 'temp';
%fldnm = 'saln';
%fldnm = 'uvel'; % there is no u/v in relax files
%fldnm = 'vvel';

nlrs = 30;     % depth layers in HYCOM TSIS

dnmb = datenum(2011,5,1);  % interpolation date
DV = datevec(dnmb);

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthrelax  = '/Net/kronos/ddmitry/hycom/TSIS/relax/';

% Read HYCOM topo:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
HH(isnan(HH))=100;
[mm,nn]=size(HH);
[JDM,IDM] = size(HH);
IJDM = IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);


%sct = 'FLYucCar';
%sct = 'Crb1';  % Alex's section
%sct = 'GoM_NAtl';
sct = 'Crb2';
switch(sct)
 case('FLYucCar');
IJ = [   434         827
         434         411
         596         371
         706         269
        1380         269];
 case('Crb1');
 IJ = [1    601
       1401 601];
 case('GoM_NAtl');
 IJ = [          48         469
         225         623
         396         662
         422         511
         535         407
         745         417
         970         236
        1345         175
        1303         849];
 case('Crb2');
  IJ = [471 290
        471 520];
end

% Section: Fl - Yucatan - Caribbean - N Atl
%IJ = [   434         827
%         434         411
%         596         371
%         706         269
%        1380         269];
%
cca = length(IJ);

SCT.End_Pnts = IJ;

IIs = [];
JJs = [];
for ii=1:cca-1
  i1 = IJ(ii,1);
  i2 = IJ(ii+1,1);
  j1 = IJ(ii,2);
  j2 = IJ(ii+1,2);
  [I,J] = sub_xsct_indx(i1,j1,i2,j2);
  I = I(:);
  J = J(:);
  if ii>1
    I=I(2:end); 
    J=J(2:end);
  end
  IIs = [IIs;I];
  JJs = [JJs;J];
end

Ind = sub2ind([mm,nn],JJs,IIs);
SCT.II  = IIs;
SCT.JJ  = JJs;
SCT.Ind = Ind;

% Compute distances along the contour
% Get bottom depth
clear dx Hb
dx(1) = 0;
ni = length(IIs);
for ip=1:ni-1
  i0 = SCT.II(ip);
  j0 = SCT.JJ(ip);
  i1 = SCT.II(ip+1);
  j1 = SCT.JJ(ip+1);
  x0 = LON(j0,i0);
  y0 = LAT(j0,i0);
  x1 = LON(j1,i1);
  y1 = LAT(j1,i1);
  dx(ip+1) = distance_spheric_coord(y1,x1,y0,x0);
  Hb(ip) = HH(j0,i0);
  SCT.Lon(ip,1) = x0;
  SCT.Lat(ip,1) = y0;
end
SCT.Lon(ip+1,1) = x1;
SCT.Lat(ip+1,1) = y1;

Hb(ni)=HH(j1,i1);
dst = cumsum(dx);

SCT.Distance_m = dst;
SCT.Hbottom    = Hb;

f_save=0;
if f_save==1;
  fsect=sprintf('%shycom_tsis_xsct%s.mat',pthoutp, sct);
  save(fsect,'SCT');
end


f_pltmap=0;
if f_pltmap==1
  figure(10); clf;
  hold on;

  contour(HH,[0 0],'k');
  contour(HH,[-6000:1000:-10],'Color',[0.8 0.8 0.8]);

  axis('equal');
  plot(SCT.II(:),SCT.JJ(:),'-');
end



% 
% Read 
switch(fldnm);
 case('temp');
  flnm = sprintf('relax_tem_m%2.2i%2.2i%4.4i',DV(3),DV(2),DV(1));
  rfld = 'temp';
  c1=3;
  c2=25;
 case('saln');
  flnm = sprintf('relax_sal_m%2.2i%2.2i%4.4i',DV(3),DV(2),DV(1));
  rfld = 'salin';
  c1=34.8;
  c2=37.2;
end

flnm  = sprintf('relax.m%2.2i%2.2i%4.4i',DV(3),DV(2),DV(1));
fina = sprintf('%s%s.a',pthrelax,flnm);
finb = sprintf('%s%s.b',pthrelax,flnm);

[F,n,m,l] = read_hycom(fina,finb,rfld);

[ZM,ZZ] = sub_zz_zm(fina,finb,HH);


clear AA zm
for ik=1:nlrs
  fprintf('Reading layer %i\n',ik);
% Section:
  
  AA(ik,:) = squeeze(F(ik,Ind));
  zm(ik,:) = squeeze(ZM(ik,Ind));
  zz(ik,:) = squeeze(ZZ(ik,Ind));  % layer interface depths
end
I=find(AA>1e20);
AA(I)=nan;



nint=400;
%CMP = create_colormap3v2(nint,c1,c2);
%cmp1 = CMP.colormap;
CMP = create_colormap8(nint,c1,c2);
cmp = CMP.colormap;
%
% Plotting section
xx = SCT.Distance_m*1e-3;
[XX,dmm] = meshgrid(xx,[1:nlrs]);

IX = [1:length(xx)];

figure(1);
set(gcf,'Position',[1456 706 1044  574]);
clf
axes('Position',[0.09 0.3 0.81 0.6]);
pcolor(XX,zm,AA); shading flat;
%pcolor(IX,zm,AA); shading flat

hold on;
colormap(cmp);
caxis([c1 c2]);

if strncmp(fldnm,'temp',4)
%contour(XX,zm,AA,[0:0.5:5],'k');
  contour(XX,zm,AA,[4 4],'r','Linewidth',2);
end

% Plot interfaces:
f_pltlrs = 1;
if f_pltlrs==1
  for ik=14:nlrs
    zlr=zz(ik,:);
    plot(xx,zlr,'-','Color',[0.4 0.4 0.4]);
  end
end


% Bottom:
xbv=[xx(1),xx,xx(end)];
zbv=[-9000,Hb,-9000];
patch(xbv,zbv,[0.7 0.7 0.7]);

xl1=0;
xl2=max(xx);
yl1=round(min(Hb))-10;
set(gca,'Color',[0 0 0],...
 'tickdir','out',...
 'xlim',[xl1 xl2],...
 'xtick',[0:500:xl2],...
 'ylim',[yl1 0],...
 'ytick',[-5000:1000:0],...
 'Fontsize',12);

DV = datevec(dnmb);
stl = sprintf('%s IAS HYCOM-TSIS RELAX %2.2i/%2.2i/%4.4i',fldnm,DV(3),DV(2),DV(1));
title(stl);
xlabel('Distance, km');

hc = colorbar;
set(hc,'Position',[0.92 0.3 0.022 0.6],...
       'TickLength',0.042,...
       'Fontsize',12);

btx = 'plot_HYCOMrelax_vsct.m';
bottom_text(btx,'pwd',1);











