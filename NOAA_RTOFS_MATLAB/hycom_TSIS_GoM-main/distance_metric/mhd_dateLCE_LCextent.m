% Use MHD time series
% find 1st "jump" in the score = bad timing of the LCE shedding 
% in forecasts
% Relation between the jump date and the LC extent 
% at time 0
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

TISL=0.17;
%TISL=0.10;
dlt0=2.5;  % threshold of discaring LC shedding events - jump in the mhd score
           % typically, LCE shedding causes jumps ~7 times
mhd0=0.4;  % threshold MHD abs. value, large deviation from true state
	   
ptht   = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig_mhd/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';


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
LCyy=[];
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

% Initial location of the LC:
  yy=LC(ik).hcst.yy;
  xx=LC(ik).hcst.xx;
  LCyy(ik,1)=max(yy);
end

CLR=[0 0.3 0.7; ...
     1 0.2 0;...
     0.3 0.9 0];

JT=JT(:);
JTprs=JTprs(:);
JTnps=JTnps(:);

nn=length(LCyy);
XX=[ones(nn,1),LCyy];
[b,bint,rr,rint,stats]=regress(JT,XX);
dfit=XX*b;

figure(1); clf;
axes('Position',[0.1 0.5 0.4 0.4]);
hold on;
plot(LCyy,JT,'.','Color',CLR(1,:),'Markersize',20);
plot(LCyy,dfit,'r--');
set(gca, 'Tickdir','out',...
	 'xlim',[24 28],...
	 'ylim',[0 95],...
	 'xtick',[24:28],...
	 'ytick',[0:20:100],...
	 'xgrid','on',...
	 'ygrid','on',...
	 'Fontsize',12);
xlabel('LC latitude');
ylabel('Days');
stl=sprintf('Days MHD<%4.2f, PIES %4.2f m, R^2=%4.2f, p=%5.3d',...
	    mhd0,TISL,stats(1),stats(3));
title(stl);

btx='mhd_dateLCE_LCextent.m';
bottom_text(btx,'pwd',1);


nn=length(LCyy);
XX=[ones(nn,1),LCyy];
[b,bint,rr,rint,stats]=regress(JTnps,XX);
dfit=XX*b;

figure(2); clf;
axes('Position',[0.1 0.5 0.4 0.4]);
hold on;
plot(LCyy,JTnps,'.','Color',CLR(1,:),'Markersize',20);
plot(LCyy,dfit,'r--');
set(gca, 'Tickdir','out',...
	 'xlim',[24 28],...
	 'ylim',[0 95],...
	 'xtick',[24:28],...
	 'ytick',[0:20:100],...
	 'xgrid','on',...
	 'ygrid','on',...
	 'Fontsize',12);
xlabel('LC latitude');
ylabel('Days');
stl=sprintf('Days MHD<%4.2f, noPIES %4.2f m, R^2=%4.2f, p=%5.3d',...
	    mhd0,TISL,stats(1),stats(3));
title(stl);
bottom_text(btx,'pwd',1);


% Persistence
[b,bint,rr,rint,stats]=regress(JTprs,XX);
dfit=XX*b;

figure(3); clf;
axes('Position',[0.1 0.5 0.4 0.4]);
hold on;
plot(LCyy,JTprs,'.','Color',CLR(1,:),'Markersize',20);
plot(LCyy,dfit,'r--');
set(gca, 'Tickdir','out',...
	 'xlim',[24 28],...
	 'ylim',[0 95],...
	 'xtick',[24:28],...
	 'ytick',[0:20:100],...
	 'xgrid','on',...
	 'ygrid','on',...
	 'Fontsize',12);
xlabel('LC latitude');
ylabel('Days');
stl=sprintf('Days MHD<%4.2f, Persist %4.2fm, R^2=%4.2f, p=%5.3d',...
	    mhd0,TISL,stats(1),stats(3));
title(stl);

bottom_text(btx,'pwd',1);

