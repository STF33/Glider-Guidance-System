% Calculate RMSE between NEMO and 
% and NEMO randomly selected to estimate the error
% of the random field for predictability
% See: Predictability and Information Theory. Part I: Measures of Predictability
%Timothy DelSole
%Print Publication: 01 Oct 2004
%DOI: https://doi.org/10.1175/1520-0469(2004)061<2425:PAITPI>2.0.CO;2
% J Atm Science
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 1; % save mat; =2 - load saved and finish missing dates

%Z0 = -10; % discard close to coastline points
Z0 = -200;
YR = 2011;  % run for 4 months in 2011 and 2012 to match f/casts

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';

% Min / max possible day
day_min = datenum(2010,1,1);
day_max = datenum(2012,12,31);
rnd_state = [];
Nrnd = 20;   % how many random days per 1 NEMO day for computing RMSE
dltD = 90;   % +/- # of days to be excluded around RMSE day to generate random days


% Date where RMSE to compute
if YR == 2011
  mo1 = 5;
  mo2 = 9;
elseif YR == 2012
  mo1 = 1;
  mo2 = 5;
end
dd1 = datenum(YR,mo1,1);
dd2 = datenum(YR,mo2,1)+31;
TM = [dd1:dd2]';
nstp = length(TM);
  

% Load NEMO 2 hycom indices:
fmatindx = sprintf('%snemo2hycom_indx.mat',pthmat);
fprintf('Loading %s\n',fmatindx);
load(fmatindx);


%Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
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
Iocn = find(INH==1 & HH<Z0);
Nocn=length(Iocn);
clear XM YM


if m~=INDX.HYCOM_size(1) | n~=INDX.HYCOM_size(2)
  fprintf('Saved domain: %i %i, this domain: %i %i\n',...
   INDX.HYCOM_size(1:2), m, n)
  error('Saved NEMO indices are for different HYCOM domain');
end


IJhycom = INDX.HYCOM_ocean;  % saved HYCOM pnts with NEMO indices

% 
[LONN,LATN,INN] = sub_getNEMO_info;
IXnemo = []; 

% Save the state of the random number generator
if isempty(rnd_state)
  rnd_state = rng;
else
  rng(rnd_state);
end

fcstname = 'NEMO_random';
frmseout = sprintf('%sRMSE_%s_%i.mat',pthout,fcstname,YR);

fprintf('\n\n %s %s\n',fcstname,frmseout);
fprintf(' NEMO Random RMSE: Time %s-%s\n',datestr(TM(1)),datestr(TM(end)));

clear RMSERR
RMSERR            = struct;
RMSERR.Name       = 'Random NEMO SSH RMSE';
RMSERR.Name_short = 'Random_NEMO';
RMSERR.Indx_hycom = Iocn;
RMSERR.HYCOM_dim  = [m,n];
RMSERR.TM         = TM;

% Continue from the last saved record
% skip finished mat files
last_rec=0;
if f_mat==2
  fprintf('Continue from the last unfinished mat file\n');
  if exist(frmseout,'file');
    fprintf('Loading %s\n',frmseout);
    load(frmseout);
    [a1,a2]=size(RMSERR.ERR_squared);
    last_rec=a2;
%        keyboard
  end
end

for irc=1:nstp
  dnmb0 = TM(irc);
  dv0 = datevec(dnmb0);

% Generate random dates not within dltD days 
  DRND = round(rand(Nrnd,1)*(day_max-day_min)+day_min);
  dltday = abs(DRND-dnmb0);
  cc=0;
  while min(dltday) < dltD,
    cc=cc+1;
    if cc>1000, error('Endless loop random days'); end;
    I=find(dltday < dltD);
    for kk=1:length(I)
      jj=I(kk);
      DRND(jj) = round(rand*(day_max-day_min)+day_min);
    end
    dltday = abs(DRND-dnmb0);
  end

  DV = datevec(DRND);
  nrc = length(DRND);



  dmean=1;
      
  if f_mat==2
    if irc<=last_rec 
      fprintf('%s Record exist skipping %i/%i/%i\n',fcstname,DV(ii,1:3));
      continue
    else
      fprintf('Continue from ii=%i %s %i/%i/%i\n\n',ii,fcstname,DV(ii,1:3));
      f_mat=1;
      last_rec=0;
    end
  end
%
% RMSE day - get NEMO field:
  ssh0 = sub_getSSH_nemo(dv0,LONN,LATN,INN,dmean);
%
% Loop over all random days for RMSE day
  for ii=1:nrc
    tic;
    yr  = DV(ii,1);
    mo  = DV(ii,2);
    dm  = DV(ii,3);
    dnmb= DRND(ii);

    sshn = sub_getSSH_nemo(DV(ii,:),LONN,LATN,INN,dmean);
    err2 = (ssh0-sshn).^2;

% Subsample to HYCOM grid: take closes to HYCOM grid point
% subsample to HYCOM deep basin GoM
    if isempty(IXnemo),
      Inemo = INDX.I_NEMO(:,1);
      Jnemo = INDX.J_NEMO(:,1);
      IJlin = sub2ind(size(sshn),Jnemo,Inemo);
      cnm   = 0;
      fprintf(' Finding NEMO indices for HYCOM subsampled region ...\n');
      for inx=1:Nocn
        ih0=Iocn(inx);
        II=find(IJhycom == ih0);
        if isempty(II); continue; end;
        cnm = cnm+1;
        IXnemo(cnm) = IJlin(II);
      end
    end

    RMSE = []; % note this is squared error
    RMSE = err2(IXnemo); 
%
% Spatial map of RMSE;
    f_chck=0;
    if f_chck==1
      ERR=zeros(m,n)*nan;
      ERR(Iocn)=RMSE;        % squared error
%      ERR(Iocn)=sqrt(RMSE);
      figure(11); clf;
      pcolor(ERR); shading flat;
      set(gca,'xlim',[10 600],'ylim',[400 800]);
      title('NEMO RMSE^2 mapped to HYCOM');

      figure(10); clf;
      pcolor(err2); shading flat;
      title('NEMO RMSE^2');

      itst=117;
      jtst=652;
      ilx = sub2ind(size(HH),jtst,itst);
% Find HYCOM index in NEMO
      iix = find(IJhycom==ilx);
      inm = Inemo(iix);
      jnm = Jnemo(iix);

    end

    RMSE=RMSE(:);
    if ii==1
      SUMERR = RMSE;
    else
      SUMERR = SUMERR+RMSE;
    end
  end   % end random days loop

  smm = SUMERR/nrc;
  smm=smm(:);
  RMSERR.ERR_squared(:,irc)=smm;  % mean squared error over random days
%
% Spatial RMSE:
  rms = sqrt(sum(smm)/length(smm));
  fprintf('1 record %6.3f min, RMSE=%8.4g %s \n',toc/60,rms,datestr(dnmb0));

  if mod(irc,15)==0 & f_mat==1
    fprintf('Saving %s\n',frmseout);
    save(frmseout,'RMSERR');
  end
%      keyboard
end

if f_mat==1
  fprintf('Saving %s\n',frmseout);
  save(frmseout,'RMSERR');
end


