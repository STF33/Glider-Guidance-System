% Analyze LC distance metric
% LC metrics are calculated in LC_hcst_fcst.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;
startup;

close all
clear

YR1=2009;
YR2=2010;
esim='PIES';
TISL=0.17;


ptht   = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx='LC_hcst_fcst.m';

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);



YR1=2009;
YR2=2010;
fmat=sprintf('%sLC_distance_hcst_fcst_%s_%3.3icm_%i-%i.mat',...
		   pthmat,esim,round(TISL*100),YR1,YR2);
fprintf('Loading %s\n',fmat);
load(fmat);




% Plot LC and time series of MHD


LMsk=HH*0;
LMsk(HH<0)=1;
nL=length(LC);
for ik=1:nL
  fprintf('Plotting %i\n',ik);
  
  TM=LC(ik).TM;
  DV=datevec(TM);
  
% hindcast  initial
  xx=LC(ik).hcst(1).xx;
  yy=LC(ik).hcst(1).yy;
% forecast
  xf=LC(ik).fcst(1).xx;
  yf=LC(ik).fcst(1).yy;

% hindcast  end
  xxE=LC(ik).hcst(end).xx;
  yyE=LC(ik).hcst(end).yy;
  % forecast
  xfE=LC(ik).fcst(end).xx;
  yfE=LC(ik).fcst(end).yy;


  figure(ik); clf;
  axes('Position',[0.1 0.35 0.8 0.6]);
  hold on

  pcolor(LON,LAT,LMsk); shading flat;
  colormap([0 0 0; 1 1 1]);
  freezeColors;
  contour(LON,LAT,HH,[-6000:500:-100],'Color',[0.7 0.7 0.7]);
%  contour(LON,LAT,HH,[0 0],'Color',[0.2 0.2 0.2]);
  axis('equal');
  set(gca,'tickdir','out',...
          'xlim',[-94 -81],...
          'ylim',[20 30]);

  CLR=[0, 0.4, 0.9;...
       0.7, 0.3, 0;...
       0, 0.8, 1;...
       1, 0.4, 0];
  clr=CLR(1,:);
  plot(xx,yy,'k-','Linewidth',4,'Color',clr);
  clr=CLR(2,:);
  plot(xf,yf,'k-','Linewidth',3,'Color',clr);
  clr=CLR(3,:);
  plot(xxE,yyE,'k-','Linewidth',4,'Color',clr);
  clr=CLR(4,:);
  plot(xfE,yfE,'k-','Linewidth',3,'Color',clr);

  stl=sprintf('MHD LC, %s, Fcst vs Hcst, %i cm, %i/%i-%i/%i',...
	      esim,TISL*100,DV(1,1),DV(1,2),DV(end,1),DV(end,2));
  title(stl);
  
% Legend  
  axes('Position',[0.85 0.8 0.1 0.12]);
  hold on
  clr=CLR(1,:);
  x1=0.03;
  x2=0.25;
  y1=1;
  plot([x1 x2],[y1 y1],'k-','Linewidth',4,'Color',clr);
  text(x2+0.05,y1,'h/cast, t0','Fontsize',12);

  clr=CLR(3,:);
  y1=y1-0.1;
  plot([x1 x2],[y1 y1],'k-','Linewidth',4,'Color',clr);
  text(x2+0.05,y1,'h/cast, tE','Fontsize',12);
  
  clr=CLR(2,:);
  y1=y1-0.1;
  plot([x1 x2],[y1 y1],'k-','Linewidth',3,'Color',clr);
  text(x2+0.05,y1,'f/cast, t0','Fontsize',12);
  
  clr=CLR(4,:);
  y1=y1-0.1;
  plot([x1 x2],[y1 y1],'k-','Linewidth',3,'Color',clr);
  text(x2+0.05,y1,'f/cast, tE','Fontsize',12);
  
  set(gca,'xlim',[0 0.9],...
	  'ylim',[0.6 1.1],...
	  'box','on',...
	  'xtick',[],...
	  'ytick',[]);

% MHD  
  MHD=LC(ik).MHD;
  MHDp=LC(ik).MHD_prst;
  TM=LC(ik).TM;
  Td=TM-TM(1);
  
  yl1=0;
  yl2=1.01*max([1.05,max(MHD)]);
  axes('Position',[0.1 0.08 0.8 0.2]);
  hold on;
  plot(Td,MHDp,'Linewidth',1.6,'Color',CLR(1,:)); % persistence
  plot(Td,MHD,'Linewidth',2,'Color',CLR(4,:));
  set(gca,'tickdir','out',...
          'xlim',[0 length(Td)],...
          'ylim',[yl1 yl2],...
	  'xgrid','on',...
	  'ygrid','on',...
	  'Fontsize',14);
  xlabel('Forecast days');
  ylabel('MHD Score');

  axes('Position',[0.1 0.2 0.1 0.06]);
  hold on
  y1=0.3;
  x1=0.1;
  x2=0.2;
  plot([x1 x2],[y1 y1],'k-','Linewidth',2,'Color',CLR(4,:));
  text(x2+0.05,y1,'f/cast','Fontsize',14);
  y1=0.15;
  plot([x1 x2],[y1 y1],'k-','Linewidth',1.6,'Color',CLR(1,:));
  text(x2+0.05,y1,'persist.','Fontsize',14);
  set(gca,'xlim',[0 0.4],...
	  'ylim',[0. 0.45],...
	  'box','on',...
	  'xtick',[],...
	  'ytick',[],...
	  'visible','off');
  

  btx='anls_LC_dist.m';
  bottom_text(btx,'pwd',1);
 
end


