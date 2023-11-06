% Analysis: MHD scores 3-mo forecasts
% initialized from OSE hindcasts for May2009- Dec 2010
% see mhd_LCLCE_OSEhcst_fcst.m  mhd_LCLCE_OSEhcst_prst.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;
startup;

close all
clear

ptht    = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat  = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthfig  = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx='anls_OSEfcst_mhd_LCLCE.m';

ixx=0;
for isim=1:2
  if isim==1
    esim = 'PIES';
  else
    esim = 'noPIES';
  end

  for YR=2009:2010
    im1=5;
    if YR==2010; im1=1; end;
    for imo=im1:12
      fmat=sprintf('%sMHD_LCLCE_OSE_%s_fcst_%4.4i%2.2i.mat',pthmat,esim,YR,imo);
      fprintf('Loading %s\n',fmat);
      A = load(fmat);

      ixx=ixx+1;
      DV=datevec(A.TM);
      MHD(ixx).mhd = A.MHD(:,1);
      MHD(ixx).Name = sprintf('Fcst %s %4.4i%2.2i',esim,YR,imo);
      MHD(ixx).Date_str = sprintf('%2.2i/%2.2i/%4.4i',DV(1,3:-1:1));
      MHD(ixx).esim = esim;
      MHD(ixx).Time = A.TM;
      
% Persistence
      fmat=sprintf('%sMHD_LCLCE_OSE_%sprst_%4.4i%2.2i.mat',pthmat,esim,YR,imo);
      fprintf('Loading %s\n',fmat);
      A = load(fmat);

      ixx=ixx+1;
      DV=datevec(A.TM);
      MHD(ixx).mhd = A.MHD(:,1);
      MHD(ixx).Name = sprintf('Persist %s %4.4i%2.2i',esim,YR,imo);
      MHD(ixx).Date_str = sprintf('%2.2i/%2.2i/%4.4i',DV(1,3:-1:1));
      MHD(ixx).esim = sprintf('Persist %s',esim);
      MHD(ixx).Time = A.TM;
    end
  end
end


% Check dates:
nll = length(MHD);
for ixx=1:nll
  d1 = MHD(ixx).Time(end)-MHD(ixx).Time(1);
  fprintf('Check time: %i  day1 %s, Ndays=%i\n',ixx,datestr(MHD(ixx).Time(1)),d1);
end


% Weekly pools
WK=[1:7:92]; % weekly pools
nwk=length(WK)-1;

% Mean MHD by forecast-hindcast/time groups
nll = length(MHD);
POOL = struct;
for igr=1:4
  for iwk=1:nwk
    POOL(igr).WK(iwk).pw=[];
  end
end

for ill=1:nll
  esim = MHD(ill).esim;
  switch(esim)
   case('PIES');
    igr=1;
   case('Persist PIES');
    igr=2;
   case('noPIES');
    igr=3;
   case('Persist noPIES');
    igr=4;
  end

  POOL(igr).Name = esim; 
  dmm=MHD(ill).mhd;
% Weekly pools
  for iwk=1:nwk
    id1=WK(iwk);
    id2=WK(iwk+1)-1;

    pw=POOL(igr).WK(iwk).pw;
    pw=[pw;dmm(id1:id2)];
    POOL(igr).WK(iwk).pw=pw;
  end

end

nwk = length(POOL(1).WK);
for igr=1:4
  for iwk=1:nwk
    dmm = POOL(igr).WK(iwk).pw;
    mdn = nanmedian(dmm);
    p25 = prctile(dmm,25);
    p75 = prctile(dmm,75);

    POOL(igr).WK(iwk).Mdn = mdn;
    POOL(igr).WK(iwk).p25 = p25;
    POOL(igr).WK(iwk).p75 = p75;
  end
end


CLR = [0,0.3,0.7; ...
       0.,0.6,1; ...
       0.7,0.4,0; ...
       1,0.6, 0];

figure(1); clf;
set(gcf,'Position',[1573, 560, 949, 751]);
axes('Position',[0.08 0.4 0.85 0.5]);
hold on;

dx0 = 0.45;
for igr=1:2:3
  clr=CLR(igr,:);
  if igr==1
    dx = -dx0;
  else
    dx = dx0;
  end
 
  for iwk=1:nwk
    mdn = POOL(igr).WK(iwk).Mdn;
    p25 = POOL(igr).WK(iwk).p25;
    p75 = POOL(igr).WK(iwk).p75;

    ixx = iwk;
    patch([ixx+dx ixx+dx ixx ixx],[0 mdn mdn 0],clr,'Edgecolor','none');
  
% Range:
    llw = p25;
    lup = p75;
    ix0=ixx+dx/2;
    plot([ix0-0.05 ix0+0.05],[llw llw],'k-','linewidth',2);
    plot([ix0-0.05 ix0+0.05],[lup lup],'k-','linewidth',2);
    plot([ix0 ix0],[llw lup],'k-','linewidth',1.6);

  end
end

% Plot persistence:
for igr=2:2:4
  clr=CLR(igr,:);
  if igr==2
    dx = -dx0;
  else
    dx = dx0;
  end

  for iwk=1:nwk
    mdn = POOL(igr).WK(iwk).Mdn;
    p25 = POOL(igr).WK(iwk).p25;
    p75 = POOL(igr).WK(iwk).p75;

    ixx=iwk;  
    plot([ixx+dx ixx],[mdn mdn],'-','Color',clr,'Linewidth',3);
    ix0=ixx+dx/2;
    plot([ixx+dx ixx],[p25 p25],'-','Color',clr,'Linewidth',1.6);
    plot([ixx+dx ixx],[p75 p75],'-','Color',clr,'Linewidth',1.6);
    plot([ix0 ix0],[p25 p75],':','Color',clr);
  end
end

set(gca,'tickdir','out',...
        'xlim',[0.5 nwk+0.6],...
        'xtick',[1:nwk],...
        'ygrid','on',...
        'Fontsize',12);

title('Median & IQR MHD (km) LC/LCE OSE 3-mo forecasts and Persistence');

axes('Position',[0.1 0.15 0.4 0.15]);
hold;
x0=0; 
dx0=0.5;
y0=0;
dy0=0.5;

clr=CLR(1,:);
patch([x0,x0,x0+dx0, x0+dx0],[y0,y0+dy0, y0+dy0, y0],clr,'Edgecolor','none');
text(x0+1.2*dx0,y0+dy0/2,'PIES');

clr=CLR(3,:);
y0=y0+2*dy0;
patch([x0,x0,x0+dx0, x0+dx0],[y0,y0+dy0, y0+dy0, y0],clr,'Edgecolor','none');
text(x0+1.2*dx0,y0+dy0/2,'noPIES');

set(gca,'xlim',[0 3],...
        'ylim',[0 3], ...
        'visible','off');

btx = 'anls_OSEfcst_mhd_LCLCE.m';
bottom_text(btx,'pwd',1);



