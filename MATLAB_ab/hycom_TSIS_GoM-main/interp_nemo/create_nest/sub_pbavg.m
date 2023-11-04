function Pbavg = sub_pbavg;
% Pbavg is needed to estimate montg1;
% Pbavg is constant in nest file
% in experiment 022
% here Pbavg is assumed constant as well

rg = 9806;
g  = 9.806;
hg = 2^100;
Tv= 72; % topography version
ssh2m  = 1/9.806; % convert HYCOM srfhgt to ssh in m
thref = 1e-3;  % reference value of specific volume (m**3/kg)

pthtopo   = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthnest022= '/Net/bimonthly/ddmitry/GOMl0.04/nest/';

ftopo = sprintf('%s/depth_GOMl0.04_%2.2i.nc',pthtopo,Tv); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
m=mm;
n=nn;

AB='a';
fina=sprintf('%sarchv.100-103A_%s.a',pthnest022,AB);
finb=sprintf('%sarchv.100-103A_%s.b',pthnest022,AB);

[F,nn,mm,ll] = read_hycom(fina,finb,'montg1');
M=squeeze(F);
M(M>1e10)=nan;

[F,n,m,l] = read_hycom(fina,finb,'srfhgt'); 
E=squeeze(F);
E(E>1e10)=nan; % g*zeta
ssh=ssh2m*E;  % zeta --> m

% Pb average, as estimated in HYCOM:
Pbavg = 1/thref*(E-M); % Pa = kg/(m*s2)

%keyboard
return