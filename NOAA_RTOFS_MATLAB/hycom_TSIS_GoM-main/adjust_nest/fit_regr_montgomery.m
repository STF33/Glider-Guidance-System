% Estimate montg1 needed for the nest files
% by regressing SSH 
% To find the regression coefficients
% I use experiment 022
% Dmitry Dukhovskoy, COAPS FSU, 2018
%
addpath /usr/people/ddmitry/codes/MyMatlab/
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_gom04/model_nest_ob;
startup

close all 
clear

ys=2009;
pthdat='/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/freerun/';
%pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthmat = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/datamat/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';

s_mat =1;  % save coefficients

esim='freerun';
rg = 9806;
gg  = 9.806;
hg = 2^100;
huge = hg;
%Tv= 72; % topography version
ssh2m  = 1/9.806; % convert HYCOM srfhgt to ssh in m
thref = 1e-3;  % reference value of specific volume (m**3/kg)
qthref = 1/thref; % reference val. of density

fmat=sprintf('%sregression_montg_ssh_IAS003.mat',pthmat);

cc=0;
for icyc=1:1
  for iyr=2009:2009
    for iday=1:3:365
      cc=cc+1;
      YRPLT(cc,1)=iyr;
      YRPLT(cc,2)=iday;
    end
  end
end

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

%[DX,DY]=sub_dx_dy(elon,alat);


nrc=length(YRPLT);
cc=0;
for kk=1:nrc
  iyr=YRPLT(kk,1);
  iday=YRPLT(kk,2);
  HR=0;
  
  pthi=pthdat; 
  fina = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthi,iyr,iday,HR);
  finb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthi,iyr,iday,HR);
  
  if ~exist(fina,'file');
    iday=iday+1;
    fina = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthi,iyr,iday,HR);
    finb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthi,iyr,iday,HR);
  end
  
  cc=cc+1;
  dj1=datenum(iyr,1,1);
  dnmb=dj1+iday-1;
  fprintf('\n Reading %s\n,',datestr(dnmb));
  
  fld='montg1';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>=0.001*huge)=nan;
  M1=squeeze(F);

  % sea surface height (ssh) in the archv. binary is
  % "srfhgt" to get it in m: ssh = srfhgt/(thref*9806), thref=1e-3
  fld='srfhgt';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>0.001*huge)=nan;
  ssh=squeeze(F)./gg;  % ssh m

  Im=find(~isnan(M1));
  Is=find(~isnan(ssh));
  di=Im-Is;
  if max(di)>0,
    error('Montg / ssh Land masks mismatch ...');
  end
  
  MNT(:,cc)=M1(Im);
  SSH(:,cc)=ssh(Im);
  

end

SSH=SSH';
MNT=MNT';

nt=size(SSH,1)
nx=size(SSH,2);
for ik=1:nx
  if mod(ik,1000)==0,
    fprintf('Regression %5.2f%% done ...\n',ik/nx*100);
  end
  
  X=ones(nt,1);
  X=[X,SSH(:,ik)];
  Y=MNT(:,ik);

%  bb = regress(Y,X);
  [bb,bint,r,rint,sts] = regress(Y,X);
  
  BB(ik,1)=bb(1);
  BB(ik,2)=bb(2);
  BB(ik,3)=sts(1);  % r2
end

RR=HH*nan;
RR(Im)=BB(:,3);

B0=HH*nan;
B0(Im)=BB(:,1);

B1=HH*nan;
B1(Im)=BB(:,2);

H=HH(Im);

% Correct OBs - too much noise, 
% due to relaxation zone
% - ssh and montg1 may not match well
% in the simulation
% Comparing several rows from the OB
% regression coefficient are very similar 
[m,n]=size(B1);
% Fix South OB
a1=B1(1,:);
a2=B1(2,:);
I1=find(~isnan(a1));
I1n=find(isnan(a1));
I2=find(~isnan(a2));
B1(1,I2)=a2(I2);
B1(1,I1n)=nan;

B0(1,I2)=B0(2,I2);
B0(1,I1n)=nan;

% Fix E:
a1=B1(:,n-1);
a2=B1(:,n-2);
I1=find(~isnan(a1));
I1n=find(isnan(a1));
I2=find(~isnan(a2));
B1(I2,n-1)=a2(I2);
B1(I1n,n-1)=nan;

B0(I2,n-1)=B0(I2,n-2);
B0(I1n,n-1)=nan;

% Fix. N;
a1=B1(m-1,:);
a2=B1(m-2,:);
I1=find(~isnan(a1));
I1n=find(isnan(a1));
I2=find(~isnan(a2));
B1(m-1,I2)=a2(I2);
B1(m-1,I1n)=nan;

B0(m-1,I2)=B0(m-2,I2);
B0(m-1,I1n)=nan;


fplt=0;
if fplt>0
  figure(1);
  pcolor(RR); shading flat;
  title('R^2 of regression: mntg=alf+beta*SSH, 1992-1995, 022')
  caxis([0.5 1]);

  figure(2); clf;
  axes('Position',[0.2 0.54 0.45 0.38]);
  pcolor(B0); shading flat;
  title('alf, regression: mntg=alf+beta*SSH, 1992-1995, 022')
  colorbar

  axes('Position',[0.2 0.05 0.45 0.38]);
  pcolor(B1); shading flat;
  title('beta, regression: mntg=alf+beta*SSH, 1992-1995, 022')
  colorbar
%  stt='prepare_OB_anomaly/fit_regr_montgomery.m';
  btx='fit_regr_montgomery';
  bottom_text(btx,'pwd',1);
end


fplt2=0;
if fplt2>0
  pthm='/Net/ocean/ddmitry/HYCOM/GoM/GOMl0.04/080/data_anls/';
  fnc=sprintf('%sregression_montg_ssh.m',pthm);
  fprintf('Loading %s\n',fnc);
  load(fnc);
  ab1=RGR.B1;
  ab0=RGR.B0;
  figure(5); clf;
  subplot(2,1,1);
  pcolor(ab0); shading flat;
  title('GOMl0.04-081 intrcpt');
  subplot(2,1,2);
  pcolor(ab1); shading flat;
  title('GOMl0.04-081 slope');
end

% Reconstruct montg1 from ssh output 
% from the model
f_chck=0;
if f_chck>0
  iyr=2009;
  iday=60;
  pthi = sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%i/',iyr);
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
  
  figure(10); clf;
  pcolor(eM1); shading flat;
  caxis([-6 6]);
  colorbar;
  title('Estimated montg1, regression');
  dd=sprintf('%i, day %3.3i',iyr,iday);
  text(50,340,dd);
  
  txtb='fit_regr_montgomery.m';
  bottom_text(txtb,'pwd',1);
  
  figure(11); clf;
  pcolor(M); shading flat;
  caxis([-6 6]);
  colorbar
  title(sprintf('montg1 %s',fina));
  dd=sprintf('%i, day %3.3i',iyr,iday);
  text(50,340,dd);
  
  txtb='fit_regr_montgomery.m';
  bottom_text(txtb,'pwd',1);
  
end


RGR.file       = 'fit_regr_montgomery.m';
RGR.year_start = YRPLT(1,1);
RGR.year_end   = YRPLT(end,1);
RGR.expt_used  = 'IAS0.03';
RGR.ssh_units  = 'm';
RGR.intercept  = B0;
RGR.slope      = B1;

if s_mat
  fprintf('Saving %s\n',fmat);
  save(fmat,'RGR');
end

