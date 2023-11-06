% Compare HYCOM data assimilative run 2011- June 2012 vs
% NEMO simulations
% Observations are created from NEMO fields
% use LC and LCE contours
%
% Calculate MHD for the LC contours
% for the HYCOM hindcasts and NEMO
% extracted in hycom_TSIS/extr_lc_hycom_nemo.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

%f_mhd = 1;  % = 1 - calculate MHD, =
IEXX = zeros(9,1);
IEXX(6:9) = 1;  % hindcast free run # expt
%ixx = 3;  % hindcast/free run # expt

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
btx = 'mhd_LCLE_OSSEhndcst_hycom_nemo.m';

% HYCOM-TSIS hindcast experiments:
fhnd = 'hycom_tsis_expts.mat';
load(fhnd);
%Nruns = length(EXPT);
Nruns = length(IEXX);

for ii=1:length(IEXX)
  if IEXX(ii)==0
    fprintf('%i : OFF    %s \n',ii,EXPT(ii).Name);
  else
    fprintf('%i : ON ---> %s \n',ii,EXPT(ii).Name);
  end
end

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;


% YC point:
x0 = -85;
y0 = 22;
%if f_mhd==1
% Array with NEMO LC
%  fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat); OLD
fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat);
fprintf('Loading %s\n',fmatout);
load(fmatout);
LCN   = LCXY;      % LC contour
LCEN  = LCE;       % LCEs
TMN   = LCN(1).TM;
Nlcec = length(LCEN);

nrc = length(LCN(1).XY);

for ixx=1:Nruns
 nmexp = EXPT(ixx).Name;
 pthd1 = EXPT(ixx).path;

 if IEXX(ixx)==0, continue; end; 

 fmatout = sprintf('%shycom_LCcontour_%2.2i.mat',pthmat,ixx);
 fprintf('Loading %s\n',fmatout);
 load(fmatout);

 TM  = LCXY.TM;

 fprintf('Calculating MHD %s\n',nmexp);
 for irc=1:nrc
  if mod(irc,100)==0,
   fprintf('  %6.2f%% done ...\n',irc/nrc*100);
  end

  dm1= TMN(irc);
  dnmb = TMN(irc);
  ihc= find(TM==dm1);

% NEMO NR:
  xn = LCN(1).XY(irc).X; % NEMO LC contour
  yn = LCN(1).XY(irc).Y;

% Hindcast:
  if ~isempty(ihc);
    xh1 = LCXY.XY(ihc).X;
    yh1 = LCXY.XY(ihc).Y;

    Nlc1 = length(LCEN);
    Nlc2 = length(LCE);

% Add/remove 1 LCE based on MHD 
% if adding LCE imporves MHD - add, otherwise not
    MHD_LCLCE = [];
    for ilce=1:Nlc1+1   % LC/LCE combination for hindcast OSE
      Neddy = ilce-1;  % Neddy=0 - only LC
      [Xc,Yc,flgC] = sub_combineLCLCE_EddyN(xn,yn,LCEN,irc,Neddy);

      for jlce=1:Nlc2+1   % adding eddies to the f/cast
        if ilce>1 & flgC==0
          MHD_LCLCE(jlce,ilce)= 1e20; % no LCEs for this time recrod
          continue;
        end
        Neddy = jlce-1;
        [Xhc,Yhc,flg] = sub_combineLCLCE_EddyN(xh1,yh1,LCE,ihc,Neddy);
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

    mhd_fcst = min(min(MHD_LCLCE));
    [j0,i0] = find(MHD_LCLCE==mhd_fcst,1);

    MHD(irc,1) = mhd_fcst;
    MHD(irc,2) = j0-1;     % how many eddies used for LC/LCE contour - HYCOM
    MHD(irc,3) = i0-1;     % # of eddies used for MHD, NEMO

% Plot
    f_chck = 0;
    if irc == 4 & f_chck==1
      figure(10); clf;
      set(gcf,'Position',[1705         682         801         639]);
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
      esim = EXPT(ixx).Name;
      stl = sprintf('%s %s fcst day=%i, MHD=%6.3f',esim,dstr,irc,mhd_fcst);
      title(stl);
      bottom_text(btx,'pwd',1);
      keyboard
    end

  else
   MHD(irc,1) = nan;
  end;

%      fprintf('mhd = %8.2f\n',mhd1);
%keyboard
 end
%keyboard    
 fprintf('  min/max MHD=%8.2f  %8.2f\n',min(min(MHD)),max(max(MHD)));
 fmat1 = sprintf('%sMHD_LCLCE_OSSEhndcst_hycom%2.2i_vs_nemo.mat',pthmat,ixx);
 fprintf('Saving %s\n',fmat1);
 save(fmat1,'MHD','TMN');
 
end;
%end

fprintf('All done\n');

keyboard


Td=TM-TM(1)+1;
Td=Td(:);
DV=datevec(TM);

cc=0;
for iyr=2011:2012
  for im=1:12
    i1=min(find(DV(:,1)==iyr & DV(:,2)==im));
    if ~isempty(i1),
      cc=cc+1;
      ttck(cc,1)=i1;
      tlbl{cc}=sprintf('%2.2i',im);
    end
  end
end


ymx = max(max(MHD));

figure(1); clf;
axes('Position',[0.08 0.45 0.85 0.47]);
hold on;
plot(MHD(:,1),'-','Color',[0 0.7 0.9],'Linewidth',2);
set(gca,'tickdir','out',...
        'xlim',[1 nrc],...
        'ylim',[0 1.1*ymx],...
        'xtick',ttck,...
        'xticklabel',tlbl,...
        'xgrid','on',...
        'ygrid','on');

title('MHD Scores HYCOM vs NEMO LC contours, 2011/2012');
%lgd = legend('HYCOM osse0','HYCOM osseE','HYCOM free','HYCOM AllSat');
%set(lgd,'position',[0.75 0.8 0.2 0.14]);
xlabel('Months');


bottom_text(btx,'pwd',1,'Position',[0.08 0.32 0.4 0.04]);





