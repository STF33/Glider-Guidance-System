% Calculate MHD
% For persistency - initial LC/LCE 
% for 90-day forecasts
% Initialized from OSE hindcasts: PIES/noPIES
%
% LC/LCE derived in LC_OSEhcst_fcst.m
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

TISL=0.17;  % target isoline for LC identification
esim='noPIES';
%esim='PIES';
%ys=2009; % 1 year at a time
YR1=2009;
YR2=2010;
%im1=5;  % forecast start month
%im2=12;
%if ys==2009, im1=5; end;
Z0 = -200;


fprintf('Persistence %4.3f m\n',TISL);

rg=9806; % convert pressure to depth, m
huge=1e20;


ptht    = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat  = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthfig  = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx='mhd_LCLCE_OSEhcst_prst.m';

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


% Control fields: OSEs 
fmatout1 = sprintf('%shycom_LCcontour_OSEhindcast_%s_2009-2011.mat',pthmat2,esim);
fprintf('Loading %s\n',fmatout1);
load(fmatout1);
LCN   = LCXY;
LCEN  = LCE;  % LCE contours
TMN   = LCN(1).TM;  %
Nlcec = length(LCEN);


LC = struct; 
nlc=0; % counter of forecasts 
for YR=YR1:YR2
  ys=YR;
  mday=[31;28;31;30;31;30;31;31;30;31;30;31];
  if mod(ys,4)==0; mday(2)=29; end;

  im1=1;
  im2=12;
  if ys==2009, im1=5; end;
 
  tic; 
  for imo=im1:im2
%    if (YR==2009 & imo~=6) | (YR>2009); continue; end;  % missed f/cast
% F/cast
    fmat=sprintf('%sLC_distance_OSE%s_fcst_%4.4i%2.2i.mat',...
        pthmat,esim,YR,imo);
    fprintf('Loading LCLCE f/cast %s\n',fmat);
    load(fmat);
    TM  = LCXY.TM;
    nrc = length(TM);

    xh1=[];
    yh1=[];
    fprintf('Calculating MHD \n');
    for irc=1:nrc       % 1 forecast 
      if mod(irc,30)==0,
        fprintf('  %6.2f%% done ...\n',irc/nrc*100);
        fprintf('  min/max MHD=%8.2f  %8.2f\n',min(MHD(:,1)),max(MHD(:,1)));
      end

      dnmb0 = TM(1); % initial date of f/cast
      dnmb = TM(irc);
      ihc = find(TMN==dnmb);

      if isempty(ihc)
        fprintf('Missing day in OSE: %s\n',datestr(dnmb));
        MHD(irc,1)=nan;
        continue;
      end

      
% OSE h/cast:
      xn = LCN(1).XY(ihc).X; % LC contours from OSE 
      yn = LCN(1).XY(ihc).Y;
% Persistence:
      if isempty(xh1)
        ihc0=ihc;
      end
      xh1 = LCN(1).XY(ihc0).X;
      yh1 = LCN(1).XY(ihc0).Y;
      LCE = LCEN;   % <--------- ????
      LCXY = LCN;

      Nlc1 = length(LCEN);
      Nlc2 = length(LCE);

%
% Add/remove 1 LCE based on MHD 
% if adding LCE imporves MHD - add, otherwise not
      MHD_LCLCE = [];
      for ilce=1:Nlc1+1   % LC/LCE combination for hindcast OSE
        Neddy = ilce-1;  % Neddy=0 - only LC
        [Xc,Yc,flgC] = sub_combineLCLCE_EddyN(xn,yn,LCEN,ihc,Neddy);

        for jlce=1:Nlc2+1   % adding eddies to the f/cast
          if ilce>1 & flgC==0
            MHD_LCLCE(jlce,ilce)= 1e20; % no LCEs for this time recrod
            continue;
          end
          Neddy = jlce-1;
          [Xhc,Yhc,flg] = sub_combineLCLCE_EddyN(xh1,yh1,LCE,ihc0,Neddy);
          if jlce>1 & flg==0 % no LCE for numbers = Neddy
            MHD_LCLCE(jlce,ilce)=1.e20;
            continue;
          end

          P = [Xc,Yc];    % HYCOMcontrol run
          Q = [Xhc,Yhc];  % HYCOM fcst LC 
          mhd = modified_hausdorff_distance(P,Q,'geo');
          MHD_LCLCE(jlce,ilce) = mhd;

        end
      end

% To reduce "jumps" when the LCE shedding timing is off
% or remote eddies get dispersed
      mhd_fcst = min(min(MHD_LCLCE));
      [j0,i0] = find(MHD_LCLCE==mhd_fcst,1);

      MHD(irc,1) = mhd_fcst;
      MHD(irc,2) = j0-1;     % how many eddies used for LC/LCE contour
      MHD(irc,3) = i0-1;

% Plot
      f_chck = 0;
      if irc==70 & f_chck==1
        figure(10); clf;
        hold on;
        plot(xn,yn,'r+');
        plot(xh1,yh1,'b.');
        for ilce=1:Nlc1
          xen = LCEN(ilce).XY(ihc).X;
          yen = LCEN(ilce).XY(ihc).Y;
          Neddy=i0-1;
          if ilce==Neddy
            plot(xen,yen,'r+');
          else
            plot(xen,yen,'k.');
          end
        end
        for ilce=1:Nlc2
          nXY = length(LCE(ilce).XY);  % how many records in LCE# ilce
          if irc>nXY; continue; end;
          xef = LCE(ilce).XY(irc).X;
          yef = LCE(ilce).XY(irc).Y;
          Neddy=j0-1;
          if ilce==Neddy
            plot(xef,yef,'b.');
          else
            plot(xef,yef,'c.');
          end
        end
        contour(LON,LAT,HH,[0 0],'k');
        axis('equal');
        set(gca,'xlim',[-98 -80],...
                'ylim',[17 31]);
        dstr = datestr(dnmb);
        stl = sprintf('%s %s Persistence fcst day=%i, MHD=%6.3f',esim,dstr,irc,mhd_fcst);
        title(stl);
        bottom_text(btx,'pwd',1);
        keyboard
      end

      dfcst=dnmb-dnmb0+1;
%      fprintf('Fcst day: %i, record=%i, mhd=%7.2f \n',dfcst,irc,mhd_fcst);

    end % loop for 1  3-mo forecast

    fprintf('1 fcst 3 month: %6.2f min\n\n',toc/60);

    if s_fmat
      fmat=sprintf('%sMHD_LCLCE_OSE_%sprst_%4.4i%2.2i.mat',...
                   pthmat,esim,YR,imo);
      fprintf('saving %s\n',fmat);
      save(fmat,'MHD','TM');
    end
  end  % loop for all months
end;  % years


  
  
  
  
  
  
  



