% compare montgomery fields
% in the original and adjusted nest files
addpath /usr/people/ddmitry/codes/MyMatlab/
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_gom04/model_nest_ob;
startup

close all 
clear

iyr=2009;
iday=5; % note every 4th day only


pthdat='/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/freerun/';
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthnest0 = '/Net/kronos/ddmitry/hycom/TSIS/';
pthnest022 = '/Net/kronos/ddmitry/hycom/TSIS/nest_files/';

hg     = 2^100;  % "huge" in HYCOM used for land masking ONLY!
huge   = hg;     % 0-depth values are not = huge
rg     = 9806;
gg  = 9.806;

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


% Old fields:
pthtmp   = sprintf('%snest_files/',pthnest0);
finbOLD = sprintf('%sarchv.%4.4i_%3.3i_00.b',pthtmp,iyr,iday);
finaOLD = sprintf('%sarchv.%4.4i_%3.3i_00.a',pthtmp,iyr,iday);
fld='montg1';
[F,n,m,l] = read_hycom(finaOLD,finbOLD,fld);
F(F>=0.001*huge)=nan;
MOLD=squeeze(F);

[F,n,m,l] = read_hycom(finaOLD,finbOLD,'srfhgt'); 
E=squeeze(F);
E(E>1e10)=nan;
ssh=E/gg;  % --> m


% New files:
pthnest = sprintf('%snest_adjusted/',pthnest0);
finaNEW = sprintf('%sarchv.%4.4i_%3.3i_00.a',pthnest,iyr,iday);
finbNEW = sprintf('%sarchv.%4.4i_%3.3i_00.b',pthnest,iyr,iday);
fld='montg1';
[F,n,m,l] = read_hycom(finaNEW,finbNEW,fld);
F(F>=0.001*huge)=nan;
MNEW=squeeze(F);

figure(1); clf;
pcolor(MOLD); shading flat
caxis([-5 5]);
colorbar;
title('OLD Montg');
text(-100,10,finaOLD,'Interpreter','none');

figure(2); clf;
pcolor(MNEW); shading flat
caxis([-5 5]);
colorbar;
title('NEW Montg');
text(-100,10,finaNEW,'Interpreter','none');

% Look at the OB:
npp=1;
figure(3); clf;
axes('Position',[0.1 0.65 0.85 0.25]);
hold on
plot(MOLD(end-npp,:));
plot(MNEW(end-npp,:));
legend('OLD','NEW','Location','Northwest');
title(sprintf('Montg at N OB, %i pnts indomain',npp));


axes('Position',[0.1 0.35 0.85 0.25]);
hold on
plot(MOLD(:,end-npp));
plot(MNEW(:,end-npp));

axes('Position',[0.1 0.05 0.85 0.25]);
hold on
plot(MOLD(npp+1,:));
plot(MNEW(npp+1,:));

% SSH at the OB
figure(4); clf;
axes('Position',[0.1 0.65 0.85 0.25]);
plot(ssh(end-npp,:));
title('SSH, North OB');

axes('Position',[0.1 0.35 0.85 0.25]);
plot(ssh(:,end-npp));

axes('Position',[0.1 0.05 0.85 0.25]);
plot(ssh(npp+1,:));




