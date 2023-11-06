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


% Hindcasts:
pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat);

fprintf('Loading %s\n',fmatout);
load(fmatout);

TM1 = LCXY(1).TM;  % Nemo
TM  = LCXY(2).TM;
TM4 = LCXY(4).TM;  % free run

nrc = length(LCXY(1).XY);
nrcf= length(LCXY(4).XY);

for irc=1:nrc
  dm1= TM1(irc);
  ihc= find(TM==dm1);
  ifr= find(TM4==dm1);

  xn = LCXY(1).XY(irc).X; % NEMO LC contour
  yn = LCXY(1).XY(irc).Y;

  if ~isempty(ihc);
				xh1 = LCXY(2).XY(ihc).X;
				yh1 = LCXY(2).XY(ihc).Y;

				xh2 = LCXY(3).XY(ihc).X;
				yh2 = LCXY(3).XY(ihc).Y;

				xh4 = LCXY(5).XY(ihc).X;
				yh4 = LCXY(5).XY(ihc).Y;

				P = [xn,yn];
				Q = [xh1, yh1];
				mhd1 = modified_hausdorff_distance(P,Q,'geo');

				Q = [xh2, yh2];
				mhd2 = modified_hausdorff_distance(P,Q,'geo');
						
    Q = [xh4, yh4];
    mhd4 = modified_hausdorff_distance(P,Q,'geo');
  else
    mhd1 = nan;
    mhd2 = nan;
    mhd4 = nan;
  end;

  if ~isempty(ifr)
    xh3 = LCXY(4).XY(irc).X;
    yh3 = LCXY(4).XY(irc).Y;

    Q = [xh3, yh3];
    mhd3 = modified_hausdorff_distance(P,Q,'geo');
  else
    mhd3=nan;
  end


  MHD(irc,1) = mhd1;
  MHD(irc,2) = mhd2;
  MHD(irc,3) = mhd3;
  MHD(irc,4) = mhd4;

end

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
plot(MHD(:,2),'-','Color',[0.9 0.5 0.0],'Linewidth',2);
plot(MHD(:,3),'-','Color',[0  0.9 0.2],'Linewidth',2);
plot(MHD(:,4),'-','Color',[1 0 0.5],'Linewidth',2);
set(gca,'tickdir','out',...
        'xlim',[1 nrc],...
        'ylim',[0 1.1*ymx],...
        'xtick',ttck,...
        'xticklabel',tlbl,...
        'xgrid','on',...
        'ygrid','on');

title('MHD Scores HYCOM vs NEMO LC contours, 2011/2012');
lgd = legend('HYCOM osse0','HYCOM osseE','HYCOM free','HYCOM AllSat');
set(lgd,'position',[0.75 0.8 0.2 0.14]);
xlabel('Months');

btx = 'mhd_osse_hindcasts_hycom_nemo.m';

bottom_text(btx,'pwd',1,'Position',[0.08 0.32 0.4 0.04]);





