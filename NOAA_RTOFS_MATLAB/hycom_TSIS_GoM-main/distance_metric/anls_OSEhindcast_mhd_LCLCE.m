% Analyze MHD for OSE hindcasts

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

% Figure to plot:
f_mhd_tser = 0; % plot MHD time series 


pthmat = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
fmat = sprintf('%sMHD_LCLCE_OSEhindcast.mat',pthmat);
btx  = 'anls_OSEhindcast_mhd_LCLCE.m';

load(fmat);

DV=datevec(TMN);
yr1=DV(1,1);
yr2=DV(end,1);
cc=0;
for yr=yr1:yr2
  for imo=1:12
    I=find(DV(:,1)==yr & DV(:,2)==imo);
%    dmm = nanmean(MHD(I));
    dmm = nanmedian(MHD(I));
    fc1 = prctile(MHD(I),25);
    fc2 = prctile(MHD(I),75);
    cc=cc+1;
    mhd_mn(cc,1)=dmm;
    mhd_mn(cc,2)=fc1;
    mhd_mn(cc,3)=fc2;
  end
end

nrec = cc;

figure(1); clf;
set(gcf,'Position',[1428         706         874         511]);
hold on;

clr=[0, 0.4, 0.8];

fprintf('Plotting ...\n');
dx=0.48;
dlx=0.1;
for icc=1:nrec
  ixx=icc;
  mhd = mhd_mn(icc,1);
  patch([ixx-dx ixx-dx ixx+dx ixx+dx],[0 mhd mhd 0],clr,'Edgecolor',[0 0 0]);
  llw = mhd_mn(icc,2);
  lup = mhd_mn(icc,3);
	plot([ixx-dlx ixx+dlx],[llw llw],'k-','linewidth',1.5);
	plot([ixx-dlx ixx+dlx],[lup lup],'k-','linewidth',1.5);
	plot([ixx ixx],[llw lup],'k-','linewidth',1.5);

end

set(gca,'tickdir','out',...
        'xlim',[1-dx nrec+dx],...
        'xtick',[0:nrec]);

title('Median and IQR for MHD LC/LCE OSE hindcasts PIES/noPIES, 2009-2011');
bottom_text(btx,'pwd',1);


  




