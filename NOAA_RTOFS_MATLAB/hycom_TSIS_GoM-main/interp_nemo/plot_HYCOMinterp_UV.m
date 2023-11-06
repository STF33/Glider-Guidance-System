% Plot U/V form NEMO-GLORYS interpolated 
% into HYCOM-TSIS horizontal grid
% but still on NEMO Z-levels
% This is to check fields before vertical interpolation
% onto HYCOM vertical grid
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear
close all

dnmb = datenum(2011,5,1);  % interpolation date
zz0 = -1.5;


DV = datevec(dnmb);

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthintrp  = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthmat    = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/datamat/';


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


% Get NEMO Z-levels:
% These are mid-grid depths
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
fprintf('Loading NEMO grid %s\n',fgrd);
load(fgrd);
ZMn = ZZN;  % NEMO
clear LONN LATN 

iz0 = max(find(ZZN>=zz0));
zz0 = ZZN(iz0);


nlrs = length(ZMn);

% Derive ZZ - interface layer depths
ZZn(1,1)=0;
for ik=1:nlrs
  zm1=ZMn(ik);
  dz=abs(2*(zm1-ZZn(ik)));
  ZZn(ik+1,1)=ZZn(ik,1)-dz;
end

% Merge two fields:
fnemo_indx = sprintf('%sNEMO_TO_HYCOM_TSIS_interp_pnts01.mat',pthdata);
fprintf('Loading %s\n',fnemo_indx);

% Find NEMO bounding indices:
NMI = load(fnemo_indx);
Inm = NMI.IndxHYCOM;   % indices covered by NEMO
[JH,IH] = ind2sub(size(HH),Inm);
jbot = min(JH);
jtop = max(JH);
ilft = min(IH);
irht = max(IH);




iyr = DV(1);
iday = dnmb-datenum(iyr,1,1)+1;

fuzlv = sprintf('%suvel_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
                 pthintrp,DV(3),DV(2),DV(1));
fvzlv = sprintf('%svvel_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
                 pthintrp,DV(3),DV(2),DV(1));


% Read Z-level U and V:
fprintf('Reading %s\n',fuzlv);
fid=fopen(fuzlv,'r');
for ik=1:iz0
  fprintf('Reading layer %i\n',ik);
  dmm=fread(fid,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fid,npad,'float32','ieee-be');  % read npad   
end
fclose(fid);
F = reshape(dmm,IDM,JDM);
F = F';
U = F;

fprintf('Reading %s\n',fvzlv);
fid=fopen(fvzlv,'r');
for ik=1:iz0
  fprintf('Reading layer %i\n',ik);
  dmm=fread(fid,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fid,npad,'float32','ieee-be');  % read npad   
end
fclose(fid);
F = reshape(dmm,IDM,JDM);
F = F';
V = F;



S = sqrt(U.^2+V.^2);

c1=0.;
c2=1;
CMP = create_colormap8(200,c1,c2);
cmp = CMP.colormap;
cnt = CMP.intervals;
nint= length(cnt);

figure(1); clf;
pcolor(S); shading flat;
caxis([c1 c2]);
colormap(cmp);

hold on;
plot([ilft irht],[jbot jbot],'r-');
plot([ilft irht],[jtop jtop],'r-');
plot([ilft ilft],[jbot jtop],'r-');
plot([irht irht],[jbot jtop],'r-');


stl=sprintf('NEMO U/V on HYCOM horiz grid, zz0=%7.1fm, %s',zz0,datestr(dnmb));
title(stl);

clb = colorbar;
set(clb,'Position',[0.92 0.1 0.02 0.8]);





