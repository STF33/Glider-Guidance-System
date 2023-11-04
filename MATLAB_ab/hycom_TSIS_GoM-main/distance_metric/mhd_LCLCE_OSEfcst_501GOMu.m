% Calculate MHD
% For persistency - initial LC/LCE 
% for 90-day forecasts
% Initialized from OSE hindcasts: PIES/noPIES
%
% LC/LCE derived:
% f/cast in LC_OSEhcst_fcst.m 
% GOM reanalysis in LC_501GOMu_rnls.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;
startup;

close all
clear

% ----------------
% Flags
% ---------------
s_fmat=1;

if s_fmat == 0
  fprintf(' ========== MHD is not saved !!! ======\n');
end

esim='noPIES';
%esim='PIES';
YR1=2009;
YR2=2010;
Z0 = -200;
lonW = -89.; % LCEs west of lonW are discarded, < -98 - include all LCEs


ptht    = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat  = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthfig  = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx='mhd_LCLCE_OSEfcst_501GOMu.m';

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);


% GoM region:
GOM=[366  489
  476  531
  583  560
  576  646
  508  827
  336  848
  204  829
  64  798
  19  746
  16  662
  12  578
  25  455
  71  382
  165  356
  281  400];

[XM,YM]=meshgrid([1:n],[1:m]);
%IN = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
Iocn = find(INH==1 & HH<Z0);
Nocn=length(Iocn);
clear XM YM

IN=INH;
Imean = INH;
Imean(HH>Z0)=0;


% Control fields: 50.1GOMu reanalysis
fmatout1 = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/LCLCEcntr_501GOMu_rnls.mat';
%sprintf('%shycom_LCcontour_OSEhindcast_%s_2009-2011.mat',pthmat2,esim);
fprintf('Loading %s\n',fmatout1);
load(fmatout1);
LCN   = LCXY;
LCEN  = LCE;  % LCE contours
TMN   = LCN(1).TM;  %
Nlcec = length(LCEN);


LC = struct; 
nlc=0; % counter of forecasts 
for YR=YR1:YR2
%  ys=YR;
  mday=[31;28;31;30;31;30;31;31;30;31;30;31];
  if mod(YR,4)==0; mday(2)=29; end;

  im1=1;
  im2=12;
  if YR==2009, im1=5; end;
 
  tic; 
% Start forecasts initialized every month:
  for imo=im1:im2
%  if (YR==2009 & imo~=6) | (YR>2009); continue; end;  % missed f/cast
% F/cast
    fmat=sprintf('%sLC_distance_OSE%s_fcst_%4.4i%2.2i.mat',...
        pthmat,esim,YR,imo);
    fprintf('Loading LCLCE f/cast %s\n',fmat);
    load(fmat);
    TM  = LCXY.TM;
    nrc = length(TM);
% Persistence contours: initial state NEMO
    xPs     = [];
    yPs     = [];
    xhPs    = [];
    yhPs    = [];
    MHD     = [];
    MHDprst = [];
    ihc0    = -100;


    xh1=[];
    yh1=[];
    fprintf('Calculating MHD \n');
    for ihc=1:nrc       % hycom forecast records
      if mod(ihc,10)==0,
        fprintf('  %6.2f%% done ...\n',ihc/nrc*100);
        fprintf('  min/max MHD=%8.2f %8.2f, min.max prst MHD=%8.2f %8.2f\n',...
         min(MHD(:,1)),max(MHD(:,1)),min(MHDprst),max(MHDprst));
      end

      dnmb0 = TM(1); % initial date of f/cast
      dnmb = TM(ihc);
      irc = find(TMN==dnmb); % corresponding time of GOMu reanalysis

      if isempty(irc)
        fprintf('Missing day in GOMu reanalysis: %s\n',datestr(dnmb));
        MHD(irc,1)=nan;
        continue;
      end
      
% GOMu reanalysis:
      xn = LCN.XY(irc).X; % LC contours from 50.1GOMu
      yn = LCN.XY(irc).Y;
% Persistence:
      if ihc0<=0
        ihc0=ihc;
      end
      xPs = LCXY.XY(ihc0).X;
      yPs = LCXY.XY(ihc0).Y;
      xh1 = LCXY.XY(ihc).X;
      yh1 = LCXY.XY(ihc).Y;

      Nlc1 = sub_numbLCE(LCEN,irc);
      Nlc2 = sub_numbLCE(LCE,ihc); % # of LCE on this date 
      [EComb, icomb] = sub_lce_comb(Nlc2);
      if Nlc2 ~= length(EComb)-1;
        error(' Check EComb - does not match Nlc2-1=%i',Nlc2-1);
      end

% 50.1 GOMu - all LCEs
% Add/remove 1 LCE based on MHD 
% if adding LCE imporves MHD - add, otherwise not
      MHD_LCLCE = [];
      flg_lce1 = [];
      for ilce=1:Nlc1+1   % LC/LCE combination for GOMu
        Neddy = ilce-1;  % Neddy=0 - only LC
        [Xc,Yc,flgC] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc,Neddy,'lonW',lonW);
        if ilce>1
          flg_lce1(ilce-1)=flgC;
        end
        xn = Xc;
        yn = Yc;
      end

% For HYCOM forecast - try all possible combinations of LC/LCEs and find the best
  % that matches the GOMu
      for jlce=1:icomb   % adding eddies to the f/cast
        Neddy = jlce-1;
  %      [Xhc,Yhc,flg] = sub_combineLCLCE_EddyN(xh1,yh1,LCE,ihc,Neddy);
        [Xhc,Yhc,flg] = sub_combineLCLCE_EComb(xh1,yh1,LCE,ihc,EComb,jlce);
        if jlce>1 & flg==0 % no LCE for numbers = Neddy
          MHD_LCLCE(jlce,1)=1.e20;
          continue;
        end

        P = [Xc,Yc];    % HYCOMcontrol run
        Q = [Xhc,Yhc];  % HYCOM fcst LC 
        mhd = modified_hausdorff_distance(P,Q,'geo');
        MHD_LCLCE(jlce,1) = mhd;
      end

      mhd_fcst = min(min(MHD_LCLCE));
      [j0,i0] = find(MHD_LCLCE==mhd_fcst,1);
      MHD_LCLCE0 = MHD_LCLCE;
      mhd_fcst0 = mhd_fcst;
      MHD(ihc,1) = mhd_fcst;

%  Persistence contour
% Add/remove LCE for persistence:
      MHD_LCLCE = [];
      NlcP = sub_numbLCE(LCE,ihc0);
      [ECombP, icombP] = sub_lce_comb(NlcP);
      for jlce=1:icombP   % adding eddies to the f/cast
  %      [Xprs,Yprs,flg] = sub_combineLCLCE_EddyN(xPs,yPs,LCE,ihc0,Neddy);
        [Xprs,Yprs,flg] = sub_combineLCLCE_EComb(xPs,yPs,LCE,ihc0,ECombP,jlce);
        if jlce>1 & flg==0 % no LCE for numbers = Neddy
          MHD_LCLCE(jlce,1)=1.e20;
          continue;
        end

        P = [Xc,Yc];      % HYCOMcontrol run
        Q = [Xprs,Yprs];  % HYCOM persist 
        mhd = modified_hausdorff_distance(P,Q,'geo');
        MHD_LCLCE(jlce,1) = mhd;
      end
      mhd_fcst = min(min(MHD_LCLCE));
      [j0,i0] = find(MHD_LCLCE==mhd_fcst,1);

      MHDprst(ihc,1) = mhd_fcst;

    end    % end 3-mo f/cast

%    keyboard

    MHD(:,2) = MHDprst;

    fprintf('  min/max MHD=%8.2f  %8.2f, min/max prst MHD=%8.2f %8.2f\n',...
            min(MHD(:,1)),max(MHD(:,1)),min(MHD(:,2)),max(MHD(:,2)));

    fmat1 = sprintf('%sMHD_LCLCE_501GOMu_persist_OSEfcst%s_%4.4i%2.2i.mat',...
                    pthmat,esim,YR,imo);


    if s_fmat == 1
      fprintf('Saving %s\n',fmat1);
      save(fmat1,'MHD','TM');
    else
     fprintf(' MHD is not saved !!!\n');
    end
  end  % month
end    % year

fprintf('All done\n');





