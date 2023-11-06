% Compare HYCOM data assimilative run 2011- June 2012 vs
% NEMO simulations
% Observations are created from NEMO fields
%
% Calculate MHD for the LC contours
% for the HYCOM hindcasts and NEMO
% for reference - free run HYCOM
% is used for 2011 only
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
IEXX(2:9) = 1;  % hindcast free run # expt
%ixx = 3;  % hindcast/free run # expt

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';

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


% YC point:
x0 = -85;
y0 = 22;
%if f_mhd==1
% Array with NEMO LC
%  fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat); OLD
  fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat);
  fprintf('Loading %s\n',fmatout);
  load(fmatout);
  LCN=LCXY;

  TMN = LCN(1).TM;  % Nemo
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
      ihc= find(TM==dm1);

      xn = LCN(1).XY(irc).X; % NEMO LC contour
      yn = LCN(1).XY(irc).Y;

      dLC = distance_spheric_coord(yn,xn,y0,x0);
      maxD = max(dLC);

      if ~isempty(ihc);
        xh1 = LCXY.XY(ihc).X;
        yh1 = LCXY.XY(ihc).Y;

        P = [xn,yn];
        Q = [xh1, yh1];
        mhd1 = modified_hausdorff_distance(P,Q,'geo');

      else
        mhd1 = nan;
      end;

      MHD(irc,1)  = mhd1;
      MAXD(irc,1)= maxD; % max dist to LC farthest point from the LC
%      fprintf('mhd = %8.2f\n',mhd1);
%keyboard
    end
%keyboard    
    fprintf('  min/max MHD=%8.2f  %8.2f\n',min(MHD),max(MHD));
    fmat1 = sprintf('%sMHD_hycom%2.2i_vs_nemo.mat',pthmat,ixx);
    fprintf('Saving %s\n',fmat1);
    save(fmat1,'MHD','TMN','MAXD');
   
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

btx = 'mhd_osse_hindcasts_hycom_nemo.m';

bottom_text(btx,'pwd',1,'Position',[0.08 0.32 0.4 0.04]);





