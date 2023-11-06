% Calculate MHD using LC/LCE contours
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;
startup;

close all
clear

ys = 2009;
ye = 2011;
esim1 = 'noPIES';
esim2 = 'PIES';


% YC point:
x0  = -85;
y0  =  22;
lnW = -100;

pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

%if f_mhd==1
% Array with NEMO LC
%  fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat); OLD
fmatout1 = sprintf('%shycom_LCcontour_OSEhindcast_%s_%2.2i-%2.2i.mat',...
                   pthmat,esim1,ys,ye);
fmatout2 = sprintf('%shycom_LCcontour_OSEhindcast_%s_%2.2i-%2.2i.mat',...
                   pthmat,esim2,ys,ye);

% OSE1:
fprintf('Loading %s\n',fmatout1);
load(fmatout1);
LCN   = LCXY;
LCEN  = LCE;  % LCE contours
TMN   = LCN(1).TM;  %
Nlcec = length(LCEN);

% OSE 2:
fprintf('Loading %s\n',fmatout2);
load(fmatout2);

TM  = LCXY.TM;
nrc = length(TMN);

fprintf('Calculating MHD \n');
for irc=1:nrc
  if mod(irc,30)==0,
    fprintf('  %6.2f%% done ...\n',irc/nrc*100);
    fprintf('  min/max MHD=%8.2f  %8.2f\n',min(MHD(:,1)),max(MHD(:,1)));
  end

  dm1= TMN(irc);
  ihc= find(TM==dm1);

  if isempty(ihc)
    MHD(irc,1)=nan;
    continue;
  end

  xn = LCN(1).XY(irc).X; % LC contours from OSE 1
  yn = LCN(1).XY(irc).Y;
% [Xc,Yc] = sub_combineLCLCE(xn,yn,LCEN,lnW,irc,0);

  xh1 = LCXY.XY(ihc).X;
  yh1 = LCXY.XY(ihc).Y;

%  [Xhc,Yhc] = sub_combineLCLCE(xh1,yh1,LCE,lnW,irc,0);

    Nlc1 = length(LCEN);
    Nlc2 = length(LCE);
%
% Add/remove LCE based on MHD 
% if adding LCE imporves MHD - add, otherwise not
  MHD_LCLCE = [];
  for ilce=1:Nlc1+1
    Neddy = ilce-1;
    [Xc,Yc,flgC] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc,Neddy);

    for jlce=1:Nlc2+1   % adding eddies to the f/cast
      if ilce>1 & flgC==0
        MHD_LCLCE(jlce,ilce)= 1e20; % no LCEs for this time recrod
        continue;
      end
      Neddy = jlce-1;
      [Xhc,Yhc,flg] = sub_combineLCLCE_EddyN(xh1,yh1,LCE,ihc,Neddy);

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

  if ihc>1 & abs(mhd_fcst-MHD(ihc-1,1))>20
    fprintf('Big dlt MHD %6.4f\n',abs(mhd_fcst-MHD(ihc-1,1)));
  end;

%keyboard
%  P = [Xc,Yc];
%  Q = [Xhc,Yhc];
%  mhd1 = modified_hausdorff_distance(P,Q,'geo');

  MHD(irc,1) = mhd_fcst;
  MHD(irc,2) = j0-1;     % how many eddies used for LC/LCE contour
  MHD(irc,3) = i0-1;
end

fprintf('  min/max MHD=%8.2f  %8.2f\n',min(MHD(:,1)),max(MHD(:,1)));
fmat1 = sprintf('%sMHD_LCLCE_OSEhindcast.mat',pthmat);
fprintf('Saving %s\n',fmat1);
save(fmat1,'MHD','TMN');

