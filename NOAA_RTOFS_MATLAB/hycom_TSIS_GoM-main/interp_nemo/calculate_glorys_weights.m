% Find surrounding GLORYS nodes
% for interpolation onto HYCOM
% Compute basis functions (polynomial coefficients)
% for bilinear interpolation 
% Saved for the whole HYCOM-TSIS domain including NEMO
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 1; % save mat; =2 - load saved and finish missing dates
s_par = 0; % parallel session
% 

dnmb = datenum(2011,5,1);  % interpolation date
DV = datevec(dnmb);


% 
pthnemo = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthdata = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthtopo = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';

% Read HYCOM topo:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);

[XM,YM]=meshgrid([1:n],[1:m]);

% Get GLORYS field:
[fglnm,flglrs] = sub_find_GLORYS_file(pthglorys,dnmb);

LONN = nc_varget(flglrs,'longitude');
LATN = nc_varget(flglrs,'latitude');
[LNN,LTN] = meshgrid(LONN,LATN);

ZZN  = nc_varget(flglrs,'depth');

%LONN = double(LONN);
%LATN = double(LATN);

%ssh = squeeze(nc_varget(flglrs,'zos'));

%keyboard

mm = size(LATN,1);
nn = size(LONN,1);
%[DXN,DYN]=sub_dx_dy(LONN,LATN);

%
% Find points outside NEMO 1km domain:
% NEMO interpolated weights/indices:
%fnemo_indx = sprintf('%sNEMO_TO_HYCOM_TSIS_interp_pnts01.mat',pthdata);
%fprintf('Loading %s\n',fnemo_indx);
%NMI = load(fnemo_indx);
%Inm = NMI.IndxHYCOM;   % indices covered by NEMO
Lmsk       = HH*0;
Lmsk(HH<0) = 1;
%Lmsk(Inm)  = 0;


Iocn = find(Lmsk>0);
IndxHYCOM = Iocn;

%keyboard

%
% Find GLRS points for all HYCOM domain:
if s_par>0
  delete(gcp('nocreate'))
  if exist('spool','var'),
    delete(spool);
  end

  spool = parpool('local');

end

fprintf('Finding GLRS points for HYCOM grid ...\n');
nI = length(Iocn);

IGLRS = zeros(nI,4); 
PHIGL = zeros(nI,4,4);  % basis functions for interp
%tic;
%nsum = zeros(nI,1);
%parfor ii=1:nI
for ii=1:nI
  atot=ii;
  if mod(atot,1000)==0,
    fprintf(' %6.2f%% processed ...\n',atot/nI*100);
  end
%  fprintf('ii=%i\n',ii);

  i0=Iocn(ii);
%  [jh,ih] = ind2sub(size(HH),i0);

  x0=LON(i0);
  y0=LAT(i0);
%
% For bilinear interpolation 
% need 4 surounding points:
  i1=max(find(LONN<=x0));
  i2=i1+1;
  j1=max(find(LATN<=y0));
  j2=j1+1;

  Igl=[i1,i2,i2,i1];
  Jgl=[j1,j1,j2,j2];
  IIgl = sub2ind([mm,nn],Jgl,Igl);
  

  f_chck=0;
  if f_chck==1
    clf;
    hold on;
    plot(x0,y0,'r*');  % HYCOM
    plot(LNN(IIgl),LTN(IIgl),'.');
  end

%
% Find basis functions for bilinear interpolation
  X=LNN(IIgl);
  Y=LTN(IIgl);
  XX=[1, 1, 1, 1];
  XX=[XX;X;Y;X.*Y];

  DD=eye(4);
  Phi = inv(XX)*DD;

% ------------
% Check: basis functions: Phi(1) = 1 for x1,y1 and 0 for other nodes
  f_chck=0;
  if f_chck==1
				for iph=1:4
						phi1 = Phi(iph,:);
						XY = X.*Y;
						dmm(iph) = sum(Phi*[1; X(iph); Y(iph); XY(iph)]); % should be ones
				end
% Interpolate f into xx0,yy0
				ff=[0.5; 0.2; 1; 5];
				xx0=mean(X);
				yy0=mean(Y);
    fi = 0;
    XY = [1; xx0; yy0; xx0*yy0];
    Fi = ff'*Phi*XY;

  end
%keyboard
% -----------
  IGLRS(ii,:)   = IIgl;
  PHIGL(ii,:,:) = Phi;  % basis functions - interp coefficients in rows
end



if exist('spool','var'),
  delete(spool);
end


if f_mat==1
  fout = sprintf('%sGLRS_TO_HYCOM_TSIS_interp_pnts.mat',pthdata)
  fprintf('Saving %s\n',fout);
  IndxHYCOM = Iocn;
  save(fout,'IGLRS','IndxHYCOM','PHIGL');
end


