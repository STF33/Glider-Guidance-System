% Plot TS vertical section form NEMO
% NEMO original grid
% Color T and contour S
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

fldnm = 'temp';
Nlrs = 75;     % depth layers in NEMO

dnmb = datenum(2011,7,9);  % interpolation date
DV = datevec(dnmb);

fprintf('Interpolating: %s %s\n\n',fldnm,datestr(dnmb));

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';

% 
% Load NEMO grid:
% Get NEMO grid:
f_get_grid=0;
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
if f_get_grid==1
  [LONN,LATN,ZZN] = sub_get_NEMO_grid(dnmb);
  LONN = double(LONN);
  LATN = double(LATN);
  ZZN = double(ZZN);
  fprintf('Saving grid %s\n',fgrd);
  save(fgrd,'LONN','LATN','ZZN');
else
  fprintf('Loading NEMO grid %s\n',fgrd);
  load(fgrd);
end

HH0 = LATN*0;
SCT = sub_GoM_xsct(LONN,LATN,HH0);

nlrs = length(ZZN);
IX = SCT.Indx;
nix = length(IX);

T=zeros(nlrs,nix)*nan;
fldnemo = 'toce';
for iz0=1:nlrs
  fprintf('depth level=%i\n',iz0);
  AA = sub_get_NEMO_TS(dnmb,fldnemo,iz0);
  T(iz0,:)=AA(IX);
end

% Get HBottom:
[a1,a2]=size(T);
clear Hb
for ii=1:a2
  i1=min(find(isnan(T(:,ii))));
  Hb(ii,1) = 0.5*(ZZN(i1-1)+ZZN(i1));
end
%
% Smooth
nd=3;
for ii=1:a2
  i1=ii-nd;
  i2=ii+nd;
  i1=max([1,i1]);
  i2=min([i2,a2]);

  dmm=mean(Hb(i1:i2));
  Hbf(ii,1)=dmm;
end


T = sub_fill_bottom_nans(T);


S=T*nan;
fldnemo = 'soce';
for iz0=1:nlrs
  fprintf('depth level=%i\n',iz0);
  AA = sub_get_NEMO_TS(dnmb,fldnemo,iz0);
  S(iz0,:)=AA(IX);
end
S = sub_fill_bottom_nans(S);


% =====================
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
% =====================

XX = cumsum(SCT.dX);
%[XX,dmm] = meshgrid(dmm,[1:nlrs+1]);



cl1=4;
cl2=31;

%yl2 = min(Hb)-50;
yl2 = -1500;
fgn=1;
sub_plot_TSxsct(fgn,cmp_sst,XX,ZZN,T,S,cl1,cl2,yl2);
%set(gca,'ylim',[yl2 0]);
ZV=[yl2;Hb;yl2];
XV=[0;XX;max(XX)];
patch(XV,ZV,[0 0 0]);


sttl = sprintf('T/S (gray=35, black=36), NEMO  %s',datestr(dnmb));
title(sttl);

btx = 'plot_TSxsct_NEMO.m';
bottom_text(btx,'pwd',1);




