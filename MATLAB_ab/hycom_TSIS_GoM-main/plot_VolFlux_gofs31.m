% Plot vol fluxes from the new HYCOM2.3-TSIS
% nested within the GLBb0.08 GOFS3.1 analysis
%
% Yucatan Vol Flux
% calculated in
% calc_volFlux_nest_GOFS31_41lrs.m
% calc_volFlux_hycom23_41lrs.m
%
% Transport saved by years
% 

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_gom04/plot_binary_output;
startup

clear all
close all

esim1='gofs31_nest'; % Nest files from GLBb0.08 GOFS3.1 
esim2='freerun_gofs31'; % 

FPLT=[1,1]; % experiments to plot

%pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthmat = '/Net/kronos/ddmitry/hycom/TSIS/datamat3/';

xsct_name = 'Yucatan';

ys=2019;
ye=2019;

cc=0;
if FPLT(1)==1;
  cc=cc+1;
  Tr=sub_combine_transp_v2(pthmat,esim1,ys,ye);
  TR(cc).Name=esim1;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;
end
if FPLT(2)==1;
  cc=cc+1;
  Tr=sub_combine_transp_v2(pthmat,esim2,ys,ye);
  TR(cc).Name=esim2;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;
end
nTR=cc;


figure(1); clf;
set(gcf,'Position',[1233         580        1305         747]);
axes('Position',[0.08 0.45 0.85 0.45]);
hold on;
cc=0;
for ik=1:nTR
  YR=TR(ik).YR;
  TT=TR(ik).Ttot*1e-6;
  plot(YR,TT);
  cc=cc+1;
  ltxt{cc}=TR(ik).Name;
  stx{ik}=sprintf('%s=%3.1f, Sv',TR(ik).Name,nanmean(TT));
end

yr1 = floor(min(TR(1).YR));
yr2 = ceil(max(TR(1).YR));

lg=legend(ltxt,'Location','SouthOutside','Interpreter','none');
set(lg,'Fontsize',13,'Position',[0.8 0.3 0.09 0.07]);
set(gca,'tickdir','out',...
	'xtick',[yr1:1/12:yr2],...
	'ytick',[0:4:35],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
ylabel('Sv');
%stx{1}=sprintf('TSIS=%3.1f Sv',nanmean(TH*1e-6));
%stx{2}=sprintf('GLBb=%3.1f Sv',mean(TG*1e-6));
title('Yucatan Transport, Sv');

btx='plot_VolFlux_gofs31.m'; 
bottom_text(btx,'pwd',1,'Position',[0.05 0.18 0.4 0.05]);

axes('Position',[0.08 0.3 0.4 0.1]);
text(0,0,stx,'Fontsize',14,'Interpreter','none');
set(gca,'Visible','off');
