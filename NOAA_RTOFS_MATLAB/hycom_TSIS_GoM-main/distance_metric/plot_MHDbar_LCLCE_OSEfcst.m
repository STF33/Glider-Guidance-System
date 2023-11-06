% Plot MHD bar for monthly  MHD scores 3-mo forecasts
% initialized from OSE hindcasts for May2009- Dec 2010
% see mhd_LCLCE_OSEhcst_fcst.m  mhd_LCLCE_OSEhcst_prst.m
%
% For weekly bar diagram see anls_OSEfcst_mhd_LCLCE.m
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



% Mean MHD by forecast-hindcast/time groups
nll = length(MHD);
POOL = struct;
PRST = struct;
for ifc=1:2
  for imo=1:3
    POOL(ifc).mo(imo).pm=[];
    PRST(ifc).mo(imo).pm=[];
  end
end

ID1=[1;31;61];
ID2=[30;60;91];
for ill=1:nll
  esim = MHD(ill).esim;
  switch(esim)
   case('PIES');
    igr=1;
    ifc=1;
   case('Persist PIES');
    igr=2;
    ifc=1;
   case('noPIES');
    igr=3;
    ifc=2;
   case('Persist noPIES');
    igr=4;
    ifc=2;
  end

  dmm = MHD(ill).mhd;
% Monthly pools
  for imo=1:3
    id1=ID1(imo);
    id2=ID2(imo);
    if igr==1 || igr==3
      pm = POOL(ifc).mo(imo).pm;
      POOL(ifc).Name = esim; 
    else
      pm = PRST(ifc).mo(imo).pm;
      PRST(ifc).Name = esim; 
    end
    pm = [pm;dmm(id1:id2)];
 
    if igr==1 || igr==3
      POOL(ifc).mo(imo).pm = pm;
    else
      PRST(ifc).mo(imo).pm = pm;
    end
  end

end

% Hindcast colors 
% match with hindcast for
% comparison
CLR = [0,0.3,0.7; ...
       0.7,0.4,0];
clrp = [0.5,0.5,0.5];  % persistence
btx='plot_MHDbar_LCLCE_OSEfcst.m';


% -------------------------
% Only monthly bars
% skip weekly
% ---------------------
for ifc=1:2
  for imo=1:3
    dmm = POOL(ifc).mo(imo).pm;
    mdn = nanmedian(dmm);
    p25 = prctile(dmm,25);
    p75 = prctile(dmm,75);

    POOL(ifc).mo(imo).Mdn = mdn;
    POOL(ifc).mo(imo).p25 = p25;
    POOL(ifc).mo(imo).p75 = p75;
  end

  for imo=1:3
    dmm = PRST(ifc).mo(imo).pm;
    mdn = nanmedian(dmm);
    p25 = prctile(dmm,25);
    p75 = prctile(dmm,75);

    PRST(ifc).mo(imo).Mdn = mdn;
    PRST(ifc).mo(imo).p25 = p25;
    PRST(ifc).mo(imo).p75 = p75;
  end
end


figure(1); clf;
set(gcf,'Position',[1573, 560, 949, 751]);
axes('Position',[0.1 0.5 0.6 0.3]);
hold on;

dx0 = 0.45;
for ifc=1:2
  clr=CLR(ifc,:);
  if ifc==1
    dx = -dx0;
  else
    dx = dx0;
  end

  for imo=1:3
    mdn = POOL(ifc).mo(imo).Mdn;
    p25 = POOL(ifc).mo(imo).p25;
    p75 = POOL(ifc).mo(imo).p75;

    ixx = imo;
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
CLRP = [0.,0.6,1; ...
       1,0.6, 0];
for ifc=1:2
  clr=CLRP(ifc,:);
  if ifc==1
    dx = -dx0;
  else
    dx = dx0;
  end

  for imo=1:3
    mdn = PRST(ifc).mo(imo).Mdn;
    p25 = PRST(ifc).mo(imo).p25;
    p75 = PRST(ifc).mo(imo).p75;

    ixx=imo;
    ix0=ixx+dx/2;
    plot([ixx+dx ixx],[mdn mdn],'-','Color',clr,'Linewidth',3);
    plot([ixx+dx ixx],[p25 p25],'-','Color',clr,'Linewidth',1.6);
    plot([ixx+dx ixx],[p75 p75],'-','Color',clr,'Linewidth',1.6);
    plot([ix0 ix0],[p25 p75],':','Color',clr);
  end
end
xlabel('Forecast Months');

set(gca,'tickdir','out',...
        'xlim',[0.5 3.5],...
        'xtick',[1:3],...
        'ygrid','on',...
        'Fontsize',12);

title('Median & IQR MHD LCLCE (km) OSE 3-mo forecasts and Persistence');

bottom_text(btx,'pwd',1);



