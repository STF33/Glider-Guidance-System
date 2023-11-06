% Plot ssh subtacted area average SSH from GOMl0.04
% Plot HYCOM Analysis simulations
% NRL runs
%
% Higher res. with SSH contours
% different colormap
% oceean flow is shhown as flow lines
% submit:  matlab -nodesktop -nosplash < plot_ssh5.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

% ----------------
% Flags
% ---------------
sfig=logical(0);
f_uv=1;   % =1 - plot u,v vectors
f_cont = 1; %=1 - consecutive numbering through all years 
f_vis = 0; % =0 - keep figures in the background

regn='GOMl0.04';
%E = 023;
E = 325;
expt=sprintf('%3.3i',E);
%year=2005;
ys=2014;
ye=2015;
%cycle=1; % cycle of 18yr sim
dday = 1; % day stepping for plotting
dhr = 24; % hr stepping for plotting
ds=1;
%ms=305;
%dJ1 = datenum(ys,1,1);
%nday0=dJ1+ds-1;

id1=1;
id2=365;

rg=9806;  % convert pressure to depth, m
huge=1e20;


drnm = '/Net/ocean/ddmitry/HYCOM/GoM/GOMl0.04_analysis/';
ptht = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat  = sprintf('%s/%s/data_anls/',drnm,expt);
pthfig = sprintf('%s/%s/frames_ssh/',drnm,expt);
pthout = pthfig;

% Read the topography:
ftopo = sprintf('%sdepth_GOMl0.04_72.nc',pthtopo);
HH  = nc_varget(ftopo,'Bathymetry');
alat = nc_varget(ftopo,'Latitude');
elon = nc_varget(ftopo,'Longitude');
[mm,nn]=size(HH);
m=mm;
n=nn;

[DX,DY]=sub_dx_dy(elon,alat);

% GoM region:
GOM=[4   348
     5     2
   223     8
   256    82
   357   114
   379   128
   438   127
   427   215
   397   362
     4   348];

[XM,YM]=meshgrid([1:n],[1:m]);
IN = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
clear XM YM

c1=-0.5;
c2=0.5;

nc=50;
clB=flipud(colormap_blue(nc));
clR=colormap_red(nc);
clB(nc-3,:)=[1 1 1];
clB(nc-2,:)=[1 1 1];
clB(nc-1,:)=[1 1 1];
clB(nc,:)  =[1 1 1];
clR(4,:)   =[1 1 1];
clR(3,:)   =[1 1 1];
clR(2,:)   =[1 1 1];
clR(1,:)   =[1 1 1];
cmp=[clB;clR];
cmp=smooth_colormap(cmp,5,3);
nint=length(cmp);
cnt=(c1:(c2-c1)/nint:c2);  % fake colormap to color gray nan shallow regions
%

% Prepare info for flow-lines:
CF.cf = 1;
CF.beta = 25;
CF.v_col = [0 0 0];
% Prepare locations (indices)
% where to draw flow lines
ddh=20;
[LX,LY]=meshgrid([1:ddh:nn],[1:ddh:mm]);
hhs=HH(1:ddh:mm,1:ddh:nn);
I=find(hhs>-10);
LX(I)=nan;
LY(I)=nan;
[a1,a2]=size(LX);
np=a1*a2;
LX=reshape(LX,a1*a2,1);
LY=reshape(LY,a1*a2,1);
LX=LX(~isnan(LX));
LY=LY(~isnan(LY));
CF.initial_flow_I=LX;
CF.initial_flow_J=LY;
CF.Runge_Kutta_dT=1*86400; % time stepping in RK
CF.Runge_Kutta_Nitir=3;  % # of RK itirations

cntr=0;
for year=ys:ye
  if f_cont==0
    cntr=0; % every year start frames from 0
  end
  
  pthd = sprintf('/nexsan/GOMl0.04/GOMl0.04_%s/data/%i/',...
		 expt,year);
  
  if (mod(year,4)==0 & id2>=365 ), 
    id2=366; 
  end

  if f_vis==0;
    ff=figure('Visible','off');
  else
    figure(1); clf;
  end

  if E == 325 & year == 2014,
    id1 = 091;
  else 
    id1 = 1;
  end
  
  for iday=id1:dday:id2
    sday=sprintf('%3.3i',iday);

    hr = 0;
    fin = sprintf('%sarchv.%i_%3.3i_%2.2i_3z.nc',pthd,year,iday,hr);
    ie = exist(fin,'file');
    if ~ie
      fprintf('Missing: %s\n',fin);
      continue;
    end


    cntr=cntr+1;
    dJ1=datenum(year,1,1);
    nday = dJ1+iday-1;
    dnmb = nday+hr/24;
    DV = datevec(dnmb);
  %  date_str=datestr(nday,29);
  %  disp(date_str);
    fprintf('\n iday=%i, %4.4i/%2.2i/%2.2i:%2.2ih\n',iday, DV(1:4));
    date_str=sprintf('%4.4i/%2.2i/%2.2i:%2.2ih\n',DV(1:4));

    ssh=squeeze(nc_varget(fin,'ssh'));

  %
  % Subtract anomaly:
    dmm=ssh;
    dmm(IN==0)=nan;
%    dmm(HH>-200)=nan;
    sshM=nanmean(nanmean(dmm));
    ssh=ssh-sshM;

    if ~exist('elon');
      alat=nc_varget(fin,'Latitude');
      elon=nc_varget(fin,'Longitude');
      n=length(elon);
      m=length(alat);
    end
  %keyboard
    clf
    pcolor(elon,alat,ssh); shading flat;
    caxis([c1 c2]);
    colormap(cmp);
    hold on;
    contour(elon,alat,HH,[-200 -200],'Color',[0.8 0.8 0.8]);
    contour(elon,alat,HH,[-5000:1000:-1000],'Color',[0.7 0.7 0.7]);
    stt=sprintf('SSH anls GOMl0.04-%s: %s',expt,date_str);
    title(stt,'FontSize',14);

    axis('equal');
    set(gca,'xlim',[-98 -76.6],'ylim',[18.12 31.9]);
    set(gca,'tickdir','out','Color',[0 0 0]);
    set(gcf,'Color',[1 1 1],'InvertHardCopy','off');


    if f_uv
      F = squeeze(nc_varget(fin,'u'));
      F(F>huge)=nan;
      u=squeeze(F(1,:,:));  % 1st layer U
      F = squeeze(nc_varget(fin,'v'));
      F(F>huge)=nan;
      v=squeeze(F(1,:,:));  % 1st layer U

      dd=10;
      cf=0.3;
      beta=30;
      v_col=[0.3 0.3 0.3];
      lwd=1.;
      scl=13*(alat(2)-alat(1));

      disp('drawing vectors ...');
      for i=1:dd:n
      for j=1:dd:m
	uu=u(j,i);
	vv=v(j,i);
	if (~isnan(uu) & ~isnan(vv))
	x1=elon(j,i);
	y1=alat(j,i);
	x2=x1+uu*scl;
	y2=y1+vv*scl;
	draw_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
	end
      end
      end

    end   % if f_uv

  % ==================
  % Plot ssh
  % =================
    contour(elon,alat,ssh,[0:0.1:1],'w','linewidth',1.);
    contour(elon,alat,ssh,[-1:0.1:-0.01],'w--','linewidth',1.);

    hght=[];
    lngth=[];
    mint=10;
    fsz=10;
    bxc='w';
  %posc=[0.1 0.045 0.8 0.05];
    posc=[0.91 0.15 0.75 0.05];
    mbx=mint;
    aend=0;
    [az,axc] = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);


    drawnow
    if sfig==1
  %    fnm=sprintf('goml%s-SSH%4.4i_%s',expt,year,sday);
      fnm=sprintf('ssh%s-%4.4i',expt,cntr);
      fout=[pthout,fnm];
      fprintf('Saving %s\n',fout);
      set(gcf,'InvertHardcopy','off','Color',[1 1 1]);
      print('-djpeg','-r150',fout);
    end
%  keyboard

  end;  % for iday
end; %% year

%exit



