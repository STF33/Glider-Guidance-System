% Analyze LC distance metric
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

s_fig=1;

YR1=2009;
YR2=2010;
esim='PIES';
TISL=0.10;


ptht   = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig_mhd/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx='plot_all_LC_mhd.m';

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
% no PIES
fmat=sprintf('%sLC_distance_hcst_fcst_noPIES_%3.3icm_%i-%i.mat',...
		   pthmat,round(TISL*100),YR1,YR2);
fprintf('Loading %s\n',fmat);
load(fmat);
LCnp=LC; % noPIES

% LC = PIES
fmat=sprintf('%sLC_distance_hcst_fcst_PIES_%3.3icm_%i-%i.mat',...
		   pthmat,round(TISL*100),YR1,YR2);
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
  % forecast PIES
  xfE=LC(ik).fcst(end).xx;
  yfE=LC(ik).fcst(end).yy;

  % forecast noPIES
  xfEnp=LCnp(ik).fcst(end).xx;
  yfEnp=LCnp(ik).fcst(end).yy;


  if s_fig==1
    if ik==1
      ff=figure('Position',[1064 229 1421 1113]);
    end
    figure(ff); clf;
  else
%    figure(ik); clf;
    figure('Position',[1064 229 1421 1113]); clf;
  end
  
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
          'ylim',[20 30],...
          'Fontsize',14);

  CLR=[0, 0.3, 0.8;...
       0.3, 0.9, 0;...
       0, 0.5, 1;...
       1, 0.4, 0];
  clr=CLR(1,:);
  plot(xx,yy,'k--','Linewidth',2.5,'Color',clr); % initial position
  clr=CLR(3,:);
  plot(xxE,yyE,'k-','Linewidth',2.5,'Color',clr); % hindcast end
  clr=CLR(2,:);
  plot(xfEnp,yfEnp,'k-','Linewidth',2.5,'Color',clr);
  clr=CLR(4,:);
  plot(xfE,yfE,'k-','Linewidth',2.5,'Color',clr);

  stl=sprintf('MHD LC, 3 mo Fcst vs Hcst, %i cm, %i/%i-%i/%i',...
	      TISL*100,DV(1,1),DV(1,2),DV(end,1),DV(end,2));
  title(stl);
  
% Legend  
  axes('Position',[0.74 0.8 0.16 0.12]);
  hold on
  clr=CLR(1,:);
  x1=0.03;
  x2=0.25;
  y1=1;
  plot([x1 x2],[y1 y1],'k--','Linewidth',2.5,'Color',clr);
  text(x2+0.05,y1,'h/cast, t0','Fontsize',12);

  clr=CLR(3,:);
  y1=y1-0.1;
  plot([x1 x2],[y1 y1],'k-','Linewidth',2.5,'Color',clr);
  text(x2+0.05,y1,'h/cast, tE','Fontsize',12);
  
  clr=CLR(2,:);
  y1=y1-0.1;
  plot([x1 x2],[y1 y1],'k-','Linewidth',2.5,'Color',clr);
  text(x2+0.05,y1,'f/cast noPIES, tE','Fontsize',12);
  
  clr=CLR(4,:);
  y1=y1-0.1;
  plot([x1 x2],[y1 y1],'k-','Linewidth',2.5,'Color',clr);
  text(x2+0.05,y1,'f/cast PIES, tE','Fontsize',12);
  
  set(gca,'xlim',[0 0.9],...
	  'ylim',[0.6 1.1],...
	  'box','on',...
	  'xtick',[],...
	  'ytick',[]);

% MHD  
  MHD=LC(ik).MHD;
  MHDp=LC(ik).MHD_prst; % persistence
  MHDnp=LCnp(ik).MHD;
  MHD=MHD(:);
  MHDp=MHDp(:);
  MHDnp=MHDnp(:);
  
  TM=LC(ik).TM;
  Td=TM-TM(1)+1;
  Td=Td(:);
% Add 0 to day 0  
  MHD=[0;MHD];
  MHDp=[0;MHDp];
  Td=[0;Td];
  MHDnp=[0;MHDnp];
  
  yl1=0;
  yl2=1.01*max([1.05,max(MHD)]);
  axes('Position',[0.1 0.08 0.8 0.2]);
  hold on;
  plot(Td,MHDp,'Linewidth',2,'Color',CLR(1,:)); % persistence
  plot(Td,MHD,'Linewidth',2,'Color',CLR(4,:));   % PIES
  plot(Td,MHDnp,'Linewidth',2,'Color',CLR(2,:)); % noPIES
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
  y1=0.45;
  x1=0.1;
  x2=0.2;
  plot([x1 x2],[y1 y1],'k-','Linewidth',2,'Color',CLR(4,:));
  text(x2+0.05,y1,'f/cast PIES','Fontsize',14);
  y1=0.3;
  plot([x1 x2],[y1 y1],'k-','Linewidth',2,'Color',CLR(2,:));
  text(x2+0.05,y1,'f/cast noPIES','Fontsize',14);
  y1=0.15;
  plot([x1 x2],[y1 y1],'k-','Linewidth',1.6,'Color',CLR(1,:));
  text(x2+0.05,y1,'persist.','Fontsize',14);
  set(gca,'xlim',[0 0.4],...
	  'ylim',[0. 0.45],...
	  'box','on',...
	  'xtick',[],...
	  'ytick',[],...
	  'visible','off');
  

  bottom_text(btx,'pwd',1);
%keyboard
  if s_fig==1
%    set(gcf,'Position',[1064 229 1421 1113]);
    fgnm=sprintf('%sLC_ALL_mhd-%3.3icm_%4.4i%2.2i.png',...
		 pthfig,TISL*100,DV(1,1),DV(1,2));

    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r200',fgnm);
  end
  
end


