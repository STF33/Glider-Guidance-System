% Use MHD time series
% find 1st "jump" in the score = bad timing of the LCE shedding 
% in forecasts
%
% PIES and noPIES
% 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;
startup;

close all
clear

%s_fig=1;

YR1=2009;
YR2=2010;

%TISL=0.17;
TISL=0.10;
dlt0=2.5;  % threshold of discaring LC shedding events - jump in the mhd score
           % typically, LCE shedding causes jumps ~7 times
mhd0=0.4;  % threshold MHD abs. value, large deviation from true state
	   
ptht   = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig_mhd/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx='stat_mhd.m';

% Read topography:


YR1=2009;
YR2=2010;

% no PIES
fmat=sprintf('%sLC_distance_hcst_fcst_noPIES_%3.3icm_%i-%i.mat',...
		   pthmat,round(TISL*100),YR1,YR2);
fprintf('Loading %s\n',fmat);
load(fmat);
LCnp=LC; % noPIES

fmat=sprintf('%sLC_distance_hcst_fcst_PIES_%3.3icm_%i-%i.mat',...
		   pthmat,round(TISL*100),YR1,YR2);
fprintf('Loading %s\n',fmat);
load(fmat);
nL=length(LC);

exMHD=zeros(nL,91)*nan;
exMHDp=zeros(nL,91)*nan;
exMHDnp=zeros(nL,91)*nan;

%
% Skip large jumps > n% 
for ik=1:nL
  TM=LC(ik).TM;
  Td=TM-TM(1)+1;
  Td=Td(:);
  
  MHD=LC(ik).MHD;
  MHDp=LC(ik).MHD_prst; % persistence
  MHDnp=LCnp(ik).MHD;   % noPIES
  MHD=MHD(:);
  MHDp=MHDp(:);
  MHDnp=MHDnp(:);
%  
% Add 0 to initial state  
%  MHD=[0;MHD];
%  MHDp=[0;MHDp];
%  Td=[0;Td];
  nn=length(MHD);

  jmp=sub_mhd_jump(MHD,dlt0,mhd0);
  JT(ik)=jmp;
  jmp=sub_mhd_jump(MHDp,dlt0,mhd0);
  JTprs(ik)=jmp;
  jmp=sub_mhd_jump(MHDnp,dlt0,mhd0);
  JTnps(ik)=jmp;
  
end

[mm,nn]=size(JT);
dlt=0.2;
X1=[1:nn]-dlt;
X2=[1:nn];
X3=[1:nn]+dlt;

JT=JT(:);
JTprs=JTprs(:);
JTnps=JTnps(:);

A=[JTprs,JT,JTnps];

CLR=[0 0.3 0.7; ...
     1 0.2 0;...
     0.3 0.9 0];

close all;
figure('Position',[940 643 1561 681]); clf;
%figure(1); clf;
axes('Position',[0.08 0.5 0.85 0.4]);
hold on;
hb=bar(X1,JTprs,dlt);
set(hb,'Facecolor',CLR(1,:));  % PIES
hb=bar(X2,JT,dlt);
set(hb,'Facecolor',CLR(2,:));  % PIES
hb=bar(X3,JTnps,dlt);
set(hb,'Facecolor',CLR(3,:));  % noPIES

set(gca,'tickdir','out',...
	'xlim',[0.5 nn+0.5],...
	'ylim',[0 91],...
	'xtick',[1:nn],...
	'ytick',[0:20:91],...
	'Fontsize',14);
xlabel('Forecasts');
ylabel('Days MHD<MHD_0');
stl=sprintf('Consecutive days from t=0, MHD<%4.1f, %4.2f m',mhd0,TISL);
title(stl); 

axes('Position',[0.2 0.08 0.2 0.3]);
hbx=boxplot(A);
set(hbx,'LineWidth',2);
set(gca,'tickdir','out',...
	'ylim',[0 98],...
	'ytick',[0:20:91],...
	'xtick',[1:3],...
	'xticklabel',{'Prst','PIES','noPIES'},...
	'Fontsize',12);
title('Days MHD<MHD_0');

btx='mhd_dateLCE.m';
bottom_text(btx,'pwd',1,'Position',[0.5 0.06 0.4 0.05]);

axes('Position',[0.6 0.2 0.2 0.2]);
XV=[0.1 0.1 0.2 0.2];
YV=[0.45 0.55 0.55 0.45];
patch(XV,YV,CLR(1,:));
text(max(XV)+0.05,mean(YV),'Persist','Fontsize',14);

YV=YV-0.2;
patch(XV,YV,CLR(2,:));
text(max(XV)+0.05,mean(YV),'PIES','Fontsize',14);
  
YV=YV-0.2;
patch(XV,YV,CLR(3,:));
text(max(XV)+0.05,mean(YV),'noPIES','Fontsize',14);
  
set(gca,'xlim',[0.05 0.6],...
	'ylim',[0. 0.6],...
	'Fontsize',16,...
	'visible','off');



