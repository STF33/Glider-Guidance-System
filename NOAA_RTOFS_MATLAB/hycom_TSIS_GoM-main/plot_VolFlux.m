% Plot hindcasts (assimilated and freeruns) 
% Yucatan Vol Flux
% calculated in
% calc_volFlux_GLBb.m - Global reanalysis
% and calc_volFlux_v3.m for HYCOM-TSIS
% Transport saved by years
%
% In order to match GLBu reanalysis Yucatan Transport
% the following adjustment have been applied:
% Barotropic transport increased by 2.5 times
% Montgomery potential is regressed using GLB reanalysis fields
% mean dynamic topography in the data/observation files (used for assimilation)
% is replaced with the mean SSH map calculated from the
% GLB reanalysis
% 

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_gom04/plot_binary_output;
startup

clear all
close all

esim1='PIES'; % adjusted nest&Montg hindcast data assimilative runs with deep PIES obs
esim2='GLBu191'; % 
esim3='noPIES'; % adjusted nest hindcast data assimilative runs without deep PIES obs
esim4='freerun_noadj'; % freerun with original nest fields
esim5='PIESnoadj'; % original hindcast data assimilative runs with deep PIES obs
esim6='noPIESnoadj'; % original hindcast data assimilative runs without deep PIES obs
esim7='GLfreerun';  % GLORYS freerun
esim8='GoMu501'; % 0.04 NCODA + GoMu reanalysis

FPLT=[1,1,1,0,0,0,0,1]; % experiments to plot

pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthmatGL = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

xsct_name = 'Yucatan';

ys=2009;
ye=2011;

cc=0;
if FPLT(1)==1;
  cc=cc+1;
  Tr=sub_combine_transp(pthmat,esim1,ys,ye);
  TR(cc).Name=esim1;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;
end
if FPLT(2)==1;
  cc=cc+1;
  Tr=sub_combine_transp(pthmat,esim2,ys,ye);
  TR(cc).Name=esim2;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;
end
if FPLT(3)==1;
  cc=cc+1;
  Tr=sub_combine_transp(pthmat,esim3,ys,ye);
  TR(cc).Name=esim3;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;
end
if FPLT(4)==1;
  cc=cc+1;
  Tr=sub_combine_transp(pthmat,esim4,ys,ye);
  TR(cc).Name=esim4;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;
end
if FPLT(5)==1;
  cc=cc+1;
  Tr=sub_combine_transp(pthmat,esim5,ys,ye);
  TR(cc).Name=esim5;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;
end
if FPLT(6)==1;
  cc=cc+1;
  Tr=sub_combine_transp(pthmat,esim6,ys,ye);
  TR(cc).Name=esim6;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;
end
if FPLT(7)==1;
  cc=cc+1;
  Tr=sub_combine_transp(pthmatGL,esim7,ys,ye);
  TR(cc).Name=esim7;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;
end
if FPLT(8)==1;
  cc=cc+1;
  Tr=sub_combine_transp(pthmat,esim8,ys,ye);
  TR(cc).Name=esim8;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;
end
nTR=cc;


figure(1); clf;
set(gcf,'Position',[1224         721        1313         608])
axes('Position',[0.08 0.45 0.85 0.45]);
hold on;
cc=0;
for ik=1:nTR
  YR=TR(ik).YR;
  TT=TR(ik).Ttot*1e-6;
  plot(YR,TT,'linewidth',2);
  cc=cc+1;
  ltxt{cc}=TR(ik).Name;
  stx{ik}=sprintf('%s=%3.1f, Sv',TR(ik).Name,nanmean(TT));
end

lg=legend(ltxt,'Location','SouthOutside','Interpreter','none');
set(lg,'Fontsize',13,'Position',[0.8 0.3 0.09 0.07]);
set(gca,'tickdir','out',...
	'xtick',[2009:1/12:2012],...
	'ytick',[0:4:40],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
ylabel('Sv');
%stx{1}=sprintf('TSIS=%3.1f Sv',nanmean(TH*1e-6));
%stx{2}=sprintf('GLBb=%3.1f Sv',mean(TG*1e-6));
title('Yucatan Transport, Sv');

btx='plot_VolFlux.m';
bottom_text(btx,'pwd',1,'Position',[0.05 0.18 0.4 0.05]);

axes('Position',[0.08 0.3 0.4 0.1]);
text(0,0,stx,'Fontsize',14,'Interpreter','none');
set(gca,'Visible','off');

