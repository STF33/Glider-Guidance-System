% Use MHD time series
% Count how many days
% fcst 
% PIES and noPIES
% is better than persistence (MHD score is lower)
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

Nps=zeros(91,1);
Nnps=zeros(91,1);
Npnp=Nps;  % PIES vs no PIES
Nnpp=Nps;  % no PIES vs PIES
cntr=Nps;

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
  i1=find(MHD<MHDp);
  Nps(i1)=Nps(i1)+1;
  i1=find(MHDnp<MHDp);
  Nnps(i1)=Nnps(i1)+1;
  ll=length(MHDp);
  cntr(1:ll)=cntr(1:ll)+1;
  
% Days PIES better No PIES 
  jd=find(MHD<MHDnp);
  Npnp(jd)=Npnp(jd)+1;
  
end

% Prob:
Nps=Nps./cntr;
Nnps=Nnps./cntr;
Npnp=Npnp./cntr;
%Nnpp=Nnpp./cntr;

Nps(Nps==0)=nan;
Nnps(Nnps==0)=nan;
Npnp(Npnp==0)=nan;
%Nnpp(Nnpp==0)=nan;

CLR=[0 0.3 0.7; ...
     1 0.2 0;...
     0.3 0.9 0];


nn=length(Nnps);
dlt=0.5;
X1=[1:nn]-dlt;
X2=[1:nn]+dlt;


close all;
tic;
f1=figure('Position',[940 643 1561 681]); clf;
%figure(1); clf;
axes('Position',[0.08 0.5 0.85 0.4]);
hold on;
for jj=1:length(Nps)
  x1=jj-0.5;
  x2=jj+0.5;
  plot([x1 x2],[Nps(jj) Nps(jj)],'-','Color',CLR(2,:),'Linewidth',3);
  j0=jj-1;
  if jj>1
    x0=j0+0.5;
    y0=Nps(j0);
    y1=Nps(jj);
    if y0~=y1
      plot([x0 x1],[y0 y1],'-','Color',CLR(2,:),'Linewidth',3);
    end
  end
% No Pies  
  plot([x1 x2],[Nnps(jj) Nnps(jj)],'-','Color',CLR(3,:),'Linewidth',1.6);
  if jj>1
    x0=j0+0.5;
    y0=Nnps(j0);
    y1=Nnps(jj);
    if y0~=y1
      plot([x0 x1],[y0 y1],'-','Color',CLR(3,:),'Linewidth',1.6);
    end
  end
  
end

set(gca,'tickdir','out',...
	'xlim',[0.5 nn+0.5],...
	'ylim',[0.45 1],...
	'xtick',[1:nn],...
	'ytick',[0.1:0.1:1],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);

xlabel('Forecast Days');
ylabel('Probability');
stl=sprintf('P(MHD_{fcst}<MHD_{prst}), %4.2f m',TISL);
title(stl); 

A=[Nps,Nnps];
axes('Position',[0.2 0.08 0.2 0.3]);
hbx=boxplot(A);
set(hbx,'LineWidth',2);
set(gca,'tickdir','out',...
	'ylim',[0 1],...
	'ytick',[0:0.1:1],...
	'xtick',[1:2],...
	'xticklabel',{'PIES','noPIES'},...
	'Fontsize',12);
title('P(MHD_f<MHD_h)');

btx='mhd_fcst_vs_persist.m';
bottom_text(btx,'pwd',1,'Position',[0.5 0.06 0.4 0.05]);

figure(f1);
axes('Position',[0.6 0.2 0.2 0.2]);
hold on
plot([0.1 0.2],[0.2 0.2],'-','Color',CLR(2,:),'Linewidth',3)
text(0.25,0.2,'PIES','Fontsize',14);

plot([0.1 0.2],[0.1 0.1],'-','Color',CLR(3,:),'Linewidth',1.6)
text(0.25,0.1,'No PIES','Fontsize',14);

set(gca,'xlim',[0.05 0.6],...
	'ylim',[0. 0.4],...
	'visible','off');
toc;

% ===============
% Prob (MHD(PIES) < MHD(no PIES))
% ===============
%
f2=figure('Position',[940 643 1551 671]); clf;
%figure(1); clf;
axes('Position',[0.08 0.5 0.85 0.4]);
hold on;
hbb=bar(Npnp,0.85);
set(hbb,'Edgecolor','none','Facecolor',[0. 0.6 1]);

set(gca,'tickdir','out',...
	'xlim',[0.5 nn+0.5],...
	'ylim',[0. 1],...
	'xtick',[1:nn],...
	'ytick',[0.1:0.1:1],...
	'Fontsize',14);

xlabel('Forecast Days');
ylabel('Probability');
stl=sprintf('P(MHD_{PIES}<MHD_{noPIES}), %4.2f m',TISL);
title(stl); 

figure(f2);
axes('Position',[0.2 0.08 0.1 0.3]);
hbx=boxplot(Npnp);
set(hbx,'LineWidth',2);
set(gca,'tickdir','out',...
	'ylim',[0 1],...
	'ytick',[0:0.1:1],...
	'Fontsize',12);
title('P(MHD_{PIES})<P(MHD_{noPIES})');


bottom_text(btx,'pwd',1,'Position',[0.5 0.06 0.4 0.05]);

