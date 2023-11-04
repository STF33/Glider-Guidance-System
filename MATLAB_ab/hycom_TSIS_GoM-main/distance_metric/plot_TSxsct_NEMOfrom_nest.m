% Plot T/S vertical sections 
% across LC/LCE in the GoM 
% NEMO nature runs
% Use NEMO+GLORYS fields interpolated into HYCOM grid for specified day
% see hycom_TSIS/interp_nemo and hycom_TSIS/interp_nemo/create_nest/create_nest
% archv files will be created with NEMO/GLORYS fields
% /Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/archv.2011_190_00.a/b
% 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;
startup;

close all
clear


Modify below for NEMO


% Set flags for extracting experiments:
EXON = zeros(9,1);
EXON(2) = 1; % select expt to plot
dnmb = datenum(2011,7,9);  % date to plot

pthd1  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
pthd2  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
pthd12 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_obs/';  % all satellites included
pthd3  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/2011_GLfreerun/'; % free run
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

%fmatout = sprintf('%sLC_coord_osse_hycomV1.mat',pthmat);

btx = 'plot_TSxsct_OSSEhndcst.m';

huge = 1.e20;
rg=9806;  % convert pressure to depth, m


fhnd = 'hycom_tsis_expts.mat';
fprintf('Loading %s\n',fhnd);
load(fhnd);

Nruns = length(EXON);

for ii=1:Nruns
  if EXON(ii)==0;
    fprintf('%i : OFF    %s \n',ii,EXPT(ii).Name);
  else
    fprintf('%i : ON ---> %s \n',ii,EXPT(ii).Name);
  end
end


for ii=1:Nruns
 fprintf('%i: %s \n',ii,EXPT(ii).path);
end


% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
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


% section to plot:
SCT = sub_GoM_xsct(LON,LAT,HH);
INDs = SCT.Indx;

dv = datevec(dnmb);
yr = dv(1);
mo = dv(2);
dm = dv(3);
iday = dnmb-datenum(yr,1,1)+1;
DV = dv;


Iexpt = find(EXON==1);
for jj=1:length(Iexpt);
  ixx = Iexpt(jj);
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;
  fmatout = sprintf('%shycom_LCcontour_%2.2i.mat',pthmat,ixx);
  fprintf('%s %s\n',nmexp,fmatout);

  fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
  finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
  fin=fina;

  ie = exist(fin,'file');

  if ~ie
    fprintf('Missing: %s\n',fin);
    continue;
  end

  fprintf('Reading %s\n',fina);
  pfld = 'temp';
  [F,nn,mm,ll] = read_hycom(fina,finb,pfld);
  F(F>huge)=nan;

  T=squeeze(F(:,INDs));
  [a1,a2]=size(T);
% Prepare for plotting - add extra bogus layer
% at the bottom, Matlab won't show it
  T(end+1,1:a2)=T(end,:);

  pfld = 'salin';
  [F,nn,mm,ll] = read_hycom(fina,finb,pfld);
  F(F>huge)=nan;

  S=squeeze(F(:,INDs));
  [a1,a2]=size(S);
% Prepare for plotting - add extra bogus layer
  S(end+1,1:a2)=S(end,:);



% Prepare vertical thicknesses for plotting  
  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e10)=nan;
  F(F<0.1)=nan;
  F=F/rg;
  Dsec=squeeze(F(:,INDs));
  Dsec(Dsec==0)=nan;
% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
  clear ZZb
  Dsec(isnan(Dsec))=0;
  ZZb(1,:)=-Dsec(1,:);
  for kk=2:l
    ZZb(kk,:)=ZZb(kk-1,:)-Dsec(kk,:);
  end
% For plotting need to have #of layers +1 interfaces
% add surface layer, otherwise all values
% will be shifted 1 layer down
  [nl,npb]=size(ZZb);
  ZZ=zeros(nl+1,npb);
  ZZ(2:nl+1,:)=ZZb;
  ZZav=ZZ;

% Depths of middle of the cells:
  ZM(1,:)=0.5*ZZ(1,:);
  for kk=1:l
    ZM(kk,:)=0.5*(ZZ(kk+1,:)+ZZ(kk,:));
  end
end




