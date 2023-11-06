% Calculate MHD for the LC and LCE contours 
% and NEMO randomly selected to estimate the error
% of the random field for predictability
% See: Predictability and Information Theory. Part I: Measures of Predictability
%Timothy DelSole
%Print Publication: 01 Oct 2004
%DOI: https://doi.org/10.1175/1520-0469(2004)061<2425:PAITPI>2.0.CO;2
% J Atm Science
%
% Add/remove LCE based on best MHD
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

f_save = 1;
lonW = -89;  % LCEs west of lonW are discarded, < -98 - include all LCEs
YR   = 2011;  % run for 4 months in 2011 and 2012 to match f/casts

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';



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



% YC point:
x0  = -85;
y0  =  22;

%if f_mhd==1
% Array with NEMO LC
%  fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat); OLD
% Combine 2010 NEMO and 2011-2012 LC contours
fmatout = sprintf('%sNEMO_LCcontour2010.mat',pthmat1);
fprintf('Loading %s\n',fmatout);
load(fmatout);
LCN=LCXY;
LCEN=LCE;  % LCE contours
TMN = LCN.TM;  % Nemo
TMN=TMN(:);

fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat1);
fprintf('Loading %s\n',fmatout);
load(fmatout);
%LCXY;
%LCEN=LCE;  % LCE contours
%TMN = LCN(1).TM;  % Nemo
dmm = LCXY.TM;
dmm = dmm(:);
TMN = [TMN;dmm];
LCN.TM = TMN;
LCN.XY = [LCN.XY,LCXY.XY];

for ik=1:4
  dmm1        = LCEN(ik).TM;
  dmm2        = LCE(ik).TM;
  LCEN(ik).TM = [dmm1,dmm2];

  dmm1             = LCEN(ik).NumbLCE;
  dmm2             = LCE(ik).NumbLCE;
  LCEN(ik).NumbLCE = [dmm1,dmm2];

  dmm1        = LCEN(ik).XY;
  dmm2        = LCE(ik).XY;
  LCEN(ik).XY = [dmm1,dmm2];
end

% Save the state of the random number generator
if isempty(rnd_state)
  rnd_state = rng;
else
  rng(rnd_state);
end

clear LCXY LCE
MHD = [];

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

% Control contour for day0 from NEMO NR:
  irc0 = find(TMN==dnmb0);
% NEMO NR - all LCEs 
% Add/remove 1 LCE based on MHD 
% if adding LCE imporves MHD - add, otherwise not
% Control run = Nature Run from NEMO
  xn = LCN.XY(irc0).X; % NEMO LC contour
  yn = LCN.XY(irc0).Y;

  MHD_LCLCE = [];
  flg_lce1 = [];
  Nlc1 = sub_numbLCE(LCEN,irc0);
  for ilce=1:Nlc1+1   % LC/LCE combination for NEMO NR
    Neddy = ilce-1;  % Neddy=0 - only LC
    [Xc,Yc,flgC] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc0,Neddy,'lonW',lonW);
    if ilce>1
      flg_lce1(ilce-1)=flgC;
    end
    xn = Xc;
    yn = Yc;
  end

  if mod(irc,10)==0,
    fprintf('  %6.2f%% done ...\n',irc/nstp*100);
  end

% Loop over all random days for RMSE day
  fprintf('Calculating MHD %s\n',datestr(dnmb0));
  MHDrand = [];
  for ii=1:nrc
    tic;
    yr  = DV(ii,1);
    mo  = DV(ii,2);
    dm  = DV(ii,3);
    dnmb= DRND(ii);

%  for ihc=1:nrc
    dm1 = DRND(ii);

    ihc = find(TMN==dnmb);
% Jan 1 2010 is missing:
    if isempty(ihc);
      fprintf('Missing contour %s\n',datestr(dnmb));
      dnmb = dnmb+1;
      ihc = find(TMN==dnmb);
    end;
    xh1 = LCN.XY(ihc).X;
    yh1 = LCN.XY(ihc).Y;

    Nlc2 = sub_numbLCE(LCEN,ihc); % # of LCE on this date 
    [EComb, icomb] = sub_lce_comb(Nlc2);
    if Nlc2 ~= length(EComb)-1;
      error(' Check EComb - does not match Nlc2-1=%i',Nlc2-1);
    end

% For HYCOM forecast - try all possible combinations of LC/LCEs and find the best
% that matches the NEMO 
    for jlce=1:icomb   % adding eddies to the f/cast
      [Xhc,Yhc,flg] = sub_combineLCLCE_EComb(xh1,yh1,LCEN,ihc,EComb,jlce);
      if jlce>1 & flg==0 % no LCE 
        MHD_LCLCE(jlce,1)=1.e20;
        continue;
      end

      P = [Xc,Yc];    % NEMO NR
      Q = [Xhc,Yhc];  % HYCOM fcst LC 
      mhd = modified_hausdorff_distance(P,Q,'geo');
      MHD_LCLCE(jlce,1) = mhd;        
    end
   
    mhd_fcst = min(min(MHD_LCLCE));
    [j0,i0] = find(MHD_LCLCE==mhd_fcst,1);
    MHD_LCLCE0 = MHD_LCLCE;
    MHDrand(ii,1) = mhd_fcst;


  end    % random days loop
%keyboard 
  MHD(irc,1) = mean(MHDrand);
     
  fprintf('  min/max MHD=%8.2f  %8.2f \n',min(MHD),max(MHD));

end;

fmat1 = sprintf('%sMHD_LCLCE_nemo_random_%i.mat',pthmat,YR);
if f_save == 1 
  fprintf('Saving %s\n',fmat1);
  save(fmat1,'MHD','TM');
else
  fprintf('  NOT SAVING !!!! \n');
end
 

fprintf('All done\n');


