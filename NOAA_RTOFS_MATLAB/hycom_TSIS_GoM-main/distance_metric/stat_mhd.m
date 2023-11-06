% To generalize MHD, 
% Compute statistics of the  LC distance metric MHD
% LC metrics are calculated in LC_hcst_fcst.m
% combine both PIES and noPIES
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
%esim='PIES';
%TISL=0.17;
TISL=0.10;
%dlt0=1.8;  % threshold of discaring LC shedding events - jump in the mhd score
dlt0=2e3;  % threshold of discaring LC shedding events - jump in the mhd score
           % typically, LCE shedding causes jumps ~7 times
med=0;  % = 0 - mean score, =1 - median score
	   

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

  if dlt0<1e3  
    mhd=sub_subsample_mhd(MHD,dlt0);
  else
    mhd=MHD;
  end
  na=length(mhd);
  exMHD(ik,1:na)=mhd;

  if dlt0<1e3
    mhd=sub_subsample_mhd(MHDp,dlt0);
  else
    mhd=MHDp;
  end
  na=length(mhd);
  exMHDp(ik,1:na)=mhd;

  if dlt0<1e3
    mhd=sub_subsample_mhd(MHDnp,dlt0);
  else
    mhd=MHDnp;
  end
  na=length(mhd);
  exMHDnp(ik,1:na)=mhd;
  
end

p25=prctile(exMHD,25);
p50=prctile(exMHD,50);
p75=prctile(exMHD,75);
p25p=prctile(exMHDp,25);
p50p=prctile(exMHDp,50);
p75p=prctile(exMHDp,75);
p25np=prctile(exMHDnp,25);
p50np=prctile(exMHDnp,50);
p75np=prctile(exMHDnp,75);
%p25=min(exMHD);
%p50=prctile(exMHD,50);
%p75=max(exMHD);
%p25p=min(exMHDp);
%p50p=prctile(exMHDp,50);
%p75p=max(exMHDp);
%p25np=min(exMHDnp);
%p50np=prctile(exMHDnp,50);
%p75np=max(exMHDnp);

MN=nanmean(exMHD); 
MNp=nanmean(exMHDp);
MNnp=nanmean(exMHDnp);

% Exclude points with few data
[mm,nn]=size(exMHD);
for ik=1:nn
  aa=exMHD(:,ik);
  in=length(find(~isnan(aa)));
  if in<5, 
    p75(ik)=p50(ik);
    p25(ik)=p50(ik);
  end
  if in<4
    p50(ik)=nan;
  end
  aa=exMHDp(:,ik);
  in=length(find(~isnan(aa)));
  if in<5, 
    p75p(ik)=p50p(ik);
    p25p(ik)=p50p(ik);
  end
  if in<4
    p50p(ik)=nan;
  end
  aa=exMHDnp(:,ik);
  in=length(find(~isnan(aa)));
  if in<5, 
    p75np(ik)=p50np(ik);
    p25np(ik)=p50np(ik);
  end
  if in<4
    p50np(ik)=nan;
  end
end

xV=[1:nn,nn:-1:1];
yV=[p75,fliplr(p25)];
yVp=[p75p,fliplr(p25p)];
yVnp=[p75np,fliplr(p25np)];

CLR=[0 0.3 0.7; ...
     1 0.2 0];

%close all;
%figure('Position',[1064 229 1421 1113]); clf;
figure(1); clf;
axes('Position',[0.1 0.5 0.8 0.4]);
hold on;
%plot(nanmean(exMHDp),'-','Color',clr,'Linewidth',2);
H=fill(xV,yVp,[0.94 0.96 1]);
set(H,'edgecolor','none');

%H=fill(xV,yVnp,[0.96 1 0.96]);
%set(H,'edgecolor','none');

%H=fill(xV,yV,[1 0.96 0.96]);
%set(H,'edgecolor','none');



if med==1
  clr=CLR(1,:);
  %plot(p25p,'--','Color',clr,'Linewidth',2);
  %plot(p75p,'--','Color',clr,'Linewidth',2);
  plot(p50p,'-','Color',clr,'Linewidth',2);

  clr=CLR(2,:);
  %plot(nanmean(exMHD),'-','Color',clr,'Linewidth',2);
  plot(p50,'-','Color',clr,'Linewidth',2);
  %plot(p25,'--','Color',clr,'Linewidth',2);
  %plot(p75,'--','Color',clr,'Linewidth',2);
  plot(p50np,'-','Color',[0.3 0.9 0],'Linewidth',2);
else
  clr=CLR(1,:);
  plot(MNp,'-','Color',clr,'Linewidth',2); % persist

  clr=CLR(2,:);
  plot(MN,'-','Color',clr,'Linewidth',2);

  plot(MNnp,'-','Color',[0.3 0.9 0],'Linewidth',2);
  
end


set(gca,'tickdir','out',...
	'xlim',[0 nn],...
	'ylim',[0 1],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
xlabel('Forecast days');
ylabel('MHD Score');

if med==1
  stl=sprintf('MHD median, no LC shedding, %3.2f m cntr',TISL);
  if dlt0>=1e3
    stl=sprintf('MHD median, with LC shedding, %3.2f m cntr',TISL);
  end
else
  stl=sprintf('MHD mean, no LC shedding, %3.2f m cntr',TISL);
  if dlt0>=1e3
    stl=sprintf('MHD mean, with LC shedding, %3.2f m cntr',TISL);
  end
end


title(stl); 

  axes('Position',[0.1 0.8 0.1 0.06]);
  hold on
  y1=0.48;
  x1=0.1;
  x2=0.2;
  plot([x1 x2],[y1 y1],'k-','Linewidth',2,'Color',CLR(2,:));
  text(x2+0.05,y1,'f/cast PIES','Fontsize',14);
  y1=0.3;
  plot([x1 x2],[y1 y1],'k-','Linewidth',2,'Color',[0.3 0.9 0]);
  text(x2+0.05,y1,'f/cast noPIES','Fontsize',14);
  y1=0.12;
  plot([x1 x2],[y1 y1],'k-','Linewidth',2,'Color',CLR(1,:));
  text(x2+0.05,y1,'persist.','Fontsize',14);
  set(gca,'xlim',[0 0.4],...
          'ylim',[0.1 0.49],...
          'box','on',...
          'xtick',[],...
          'ytick',[],...
          'visible','off');


bottom_text(btx,'pwd',1,'Position',[0.04 0.38 0.4 0.1]);
  



