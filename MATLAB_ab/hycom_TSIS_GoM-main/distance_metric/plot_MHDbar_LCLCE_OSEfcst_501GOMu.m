% Plot MHD bar for monthly  MHD scores 3-mo forecasts
% initialized from OSE hindcasts for May2009- Dec 2010
% f/casts are compared against 50.1 GOMu reanalysis
% MHD computed in mhd_LCLCE_OSEfcst_501GOMu.m
% only monhtly mean statistics is plotted
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
btx     = 'plot_MHDbar_LCLCE_OSEfcst_501GOMu.m';

MHD = sub_combine_MHD_OSEfcst;

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
    lpm = length(pm);
    ldm = length(dmm);
    if ldm < id2
      dmm(ldm+1:id2) = nan;
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


% -------------------------
% Only monthly bars
% skip weekly
% ---------------------
for ifc=1:2
  for imo=1:3
    dmm = POOL(ifc).mo(imo).pm;
    mdn = nanmedian(dmm);
%    mdn = nanmean(dmm);
    p25 = prctile(dmm,25);
    p75 = prctile(dmm,75);

    POOL(ifc).mo(imo).Mdn = mdn;
    POOL(ifc).mo(imo).p25 = p25;
    POOL(ifc).mo(imo).p75 = p75;
  end

  for imo=1:3
    dmm = PRST(ifc).mo(imo).pm;
    mdn = nanmedian(dmm);
%    mdn = nanmean(dmm);
    p25 = prctile(dmm,25);
    p75 = prctile(dmm,75);

    PRST(ifc).mo(imo).Mdn = mdn;
    PRST(ifc).mo(imo).p25 = p25;
    PRST(ifc).mo(imo).p75 = p75;
  end
end


fprintf('Plotting statistics\n');
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
% Persistence:
clrP = [0.7 0.7 0.7];
for ifc=1:2
  if ifc==1
    dx = -dx0/2;
  else
    dx = dx0*3/2;
  end

  for imo=1:3
    mdn = PRST(ifc).mo(imo).Mdn;
    p25 = PRST(ifc).mo(imo).p25;
    p75 = PRST(ifc).mo(imo).p75;
    dbr = 0.02;

    ixx=imo;
    ixM=ixx+dx/2;
% Persistence:
    plot(ixM,mdn,'.','Color',clrP,'Markersize',20);
    plot([ixM-dbr ixM+dbr],[p25 p25],'-','Color',clrP,'linewidth',1.5);
    plot([ixM-dbr ixM+dbr],[p75 p75],'-','Color',clrP,'linewidth',1.5);
    plot([ixM ixM],[p25 p75],'-','Color',clrP,'linewidth',1.5);

  end
end
xlabel('Forecast Months');

set(gca,'tickdir','out',...
        'xlim',[0.5 3.5],...
        'xtick',[1:3],...
        'ygrid','on',...
        'Fontsize',12);

title('Median & IQR MHD LCLCE(km) OSE 3-mo forecasts and Persist vs 50.1GOMu');

bottom_text(btx,'pwd',1);



% Legend:
axes('Position',[0.1,0.1,0.4,0.25]);
hold on;
x1=0;
x2=1;
y1=1;
clr=CLR(1,:);
plot([x1,x2],[y1,y1],'-','Color',clr,'linewidth',20);
text(x2+0.1,y1,'PIES','Fontsize',14);

y1=2;
clr=CLR(2,:);
plot([x1,x2],[y1,y1],'-','Color',clr,'linewidth',20);
text(x2+0.1,y1,'noPIES','Fontsize',14);

set(gca,'xlim',[0 1.5],...
        'ylim',[0.8 2.3],...
        'visible','off');











