% Plot GLORYS section similar 
% to HYCOM-TSIS interpolated field
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


fid = [];

f_write = 0;
fldnm = 'temp';
%fldnm = 'saln';
%fldnm = 'uvel';
%fldnm = 'vvel';

Nlrs = 75;     % depth layers in NEMO

dnmb = datenum(2011,5,1);  % interpolation date
DV = datevec(dnmb);

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';

% Get GLORYS field:
[fglnm,flglrs] = sub_find_GLORYS_file(pthglorys,dnmb);

LONN = nc_varget(flglrs,'longitude');
LATN = nc_varget(flglrs,'latitude');
[LNN,LTN] = meshgrid(LONN,LATN);

ZZN  = nc_varget(flglrs,'depth');
ZZN = -ZZN;

%
% Get HYCOM transect
fsect=sprintf('%shycom_tsis_xsct.mat',pthoutp);
load(fsect);
SCTH = SCT;

clear SCT;
% Find section points on GLORYS grid
ni = length(SCTH.II);
cpp=0;
IIs=[];
JJs=[];
for ip=1:ni
  x0=SCTH.Lon(ip);
  y0=SCTH.Lat(ip);

  i0=max(find(LONN<=x0));
  j0=max(find(LATN<=y0));

  if isempty(IIs)
    cpp=cpp+1;
    IIs(cpp,1)=i0;
    JJs(cpp,1)=j0;
  else
% avoid repeatition
    I0 = find(IIs==i0 & JJs==j0);
    if ~isempty(I0); continue; end;
    
    cpp=cpp+1;
    IIs(cpp,1)=i0;
    JJs(cpp,1)=j0;
  end

end
Ind = sub2ind(size(LNN),JJs,IIs);
SCT.Ind = Ind;
SCT.II  = IIs;
SCT.JJ  = JJs;

% Compute distances along the contour
% Get bottom depth
clear dx 
ni = length(SCT.II);
dx(1) = 0;
ni = length(IIs);
for ip=1:ni-1
  i0 = SCT.II(ip);
  j0 = SCT.JJ(ip);
  i1 = SCT.II(ip+1);
  j1 = SCT.JJ(ip+1);
  x0 = LONN(i0);
  y0 = LATN(j0);
  x1 = LONN(i1);
  y1 = LATN(j1);
  dx(ip+1) = distance_spheric_coord(y1,x1,y0,x0);
  SCT.Lon(ip,1) = x0;
  SCT.Lat(ip,1) = y0;
end
SCT.Lon(ip+1,1) = x1;
SCT.Lat(ip+1,1) = y1;

dst = cumsum(dx);
SCT.Distance_m = dst;


switch(fldnm),
 case('temp');
  fldgl='thetao';
 case('saln');
  fldgl='so';
 case('uvel')
  fldgl='uo';
 case('vvel')
  fldgl='vo';
end

nlon = length(LONN);
nlat = length(LATN);
ndpth= length(ZZN);

for iz1=1:ndpth
  fprintf('Depth layer %i\n',iz1);
  F = squeeze(nc_varget(flglrs,fldgl,[0 iz1-1 0 0],[1 1 nlat nlon]));

%pcolor(F); shading flat;
%hold on;
%plot(IIs,JJs,'-');
%set(gca,'xlim',[1000 1600],...
%        'ylim',[1000 1400]);
%title('GLORYS section');
  AA(iz1,:) = F(Ind)';
end;

nint=400;
c1=3;
c2=25;
%CMP = create_colormap3v2(nint,c1,c2);
%cmp1 = CMP.colormap;
CMP = create_colormap8(nint,c1,c2);
cmp = CMP.colormap;
%
% Plotting section
XX = SCT.Distance_m*1e-3;

figure(1);
set(gcf,'Position',[1456    706    99     574]);
clf
axes('Position',[0.09 0.3 0.81 0.6]);
pcolor(XX,ZZN,AA); shading flat;
hold on;
colormap(cmp);
caxis([c1 c2]);

contour(XX,ZZN,AA,[0:0.5:5],'k');
contour(XX,ZZN,AA,[4 4],'k','Linewidth',2);

xl1=0;
xl2=max(XX);
yl1=-5600;
set(gca,'Color',[0.6 0.6 0.6],...
 'tickdir','out',...
 'xlim',[xl1 xl2],...
 'xtick',[0:500:xl2],...
 'ylim',[yl1 0],...
 'ytick',[-5000:1000:0],...
 'Fontsize',12);

DV = datevec(dnmb);
stl = sprintf('GLORYS %2.2i/%2.2i/%4.4i',DV(3),DV(2),DV(1));
title(stl);

hc = colorbar;
set(hc,'Position',[0.92 0.3 0.022 0.6],...
       'TickLength',0.042,...
       'Fontsize',12);

btx = 'plot_GLORYS_vsct.m';
bottom_text(btx,'pwd',1);










