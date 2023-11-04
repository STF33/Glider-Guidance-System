% Check by reconstructing montgomery pot field
% with the montg  field used to derive regression coeff.
%
addpath /usr/people/ddmitry/codes/MyMatlab/
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_gom04/model_nest_ob;
startup

close all 
clear

iyr=2009;
iday=9; % note every 4th day only

rg = 9806;
gg  = 9.806;
hg = 2^100;
huge = hg;


pthdat='/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/freerun/'; % free runs
pthi=pthdat; 
%pthi = sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%i/',iyr);
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthnest0 = '/Net/kronos/ddmitry/hycom/TSIS/';
pthnest022 = '/Net/kronos/ddmitry/hycom/TSIS/nest_files/';

hg     = 2^100;  % "huge" in HYCOM used for land masking ONLY!
huge   = hg;     % 0-depth values are not = huge
rg     = 9806;

% -------------------------
% Get grid and bath:
% My bathymetry, Model bathymetry:
% -------------------------
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
alat = nc_varget(ftopo,'mplat');
elon = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
JD=mm;
ID=nn;
m=mm;
n=nn;
HH(isnan(HH))=100;
IJDM = ID*JD;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

% Regr coeff.
fmat=sprintf('%sregression_montg_ssh_IAS003.mat',pthmat);
fprintf('Loading %s\n',fmat);
load(fmat);


fina=sprintf('%sarchv.%i_%3.3i_00.a',pthi,iyr,iday);
finb=sprintf('%sarchv.%i_%3.3i_00.b',pthi,iyr,iday);

[F,nn,mm,ll] = read_hycom(fina,finb,'montg1');
M=squeeze(F);
M(M>1e10)=nan;

[F,n,m,l] = read_hycom(fina,finb,'srfhgt'); 
E=squeeze(F);
E(E>1e10)=nan;
ssh=E/gg;  % --> m

% Estimate montg1:
B0=RGR.intercept;
B1=RGR.slope;
eM1=B0+B1.*ssh;

% Filter estimated field to get rid of 
% high-freq noise and replace the next to OB
% row with the previous one to reduce noise

nij=2; % Filter is (2*nij+1)x(2*nij+1)
eMf=sub_fltr_Gauss(nij,eM1);

% Estimated montg
figure(10); clf;
pcolor(eM1); shading flat;
caxis([-6 6]);
colorbar;
title('Reconstructed montg1, regression');
dd=sprintf('%i, day %3.3i',iyr,iday);
text(50,340,dd);

txtb='check_regr_montg.m';
bottom_text(txtb,'pwd',1);

% Actual
figure(11); clf;
pcolor(M); shading flat;
caxis([-6 6]);
colorbar
title(sprintf('montg1 %s',fina));
dd=sprintf('%i, day %3.3i',iyr,iday);
text(50,340,dd);

bottom_text(txtb,'pwd',1);

% Filtered
figure(12); clf;
pcolor(eMf); shading flat;
caxis([-6 6]);
colorbar
title(sprintf('Filtered Reconstructed montg1'));
dd=sprintf('%i, day %3.3i',iyr,iday);
text(50,340,dd);

bottom_text(txtb,'pwd',1);


npp=1;
figure(1); clf;
axes('Position',[0.1 0.65 0.85 0.25]);
hold on
plot(M(end-npp,:));
plot(eM1(end-npp,:));
plot(eMf(end-npp,:));

legend('Original','Reconstr','Filtrd');
title('Montg, North OB');

axes('Position',[0.1 0.35 0.85 0.25]);
hold on
plot(M(:,end-npp));
plot(eM1(:,end-npp));
plot(eMf(:,end-npp));

axes('Position',[0.1 0.05 0.85 0.25]);
hold on
plot(M(npp+1,:));
plot(eM1(npp+1,:));
plot(eMf(npp+1,:));


figure(2); clf;
axes('Position',[0.1 0.65 0.85 0.25]);
plot(ssh(end-npp,:));
title('SSH, North OB');

axes('Position',[0.1 0.35 0.85 0.25]);
plot(ssh(:,end-npp));

axes('Position',[0.1 0.05 0.85 0.25]);
plot(ssh(npp+1,:));




