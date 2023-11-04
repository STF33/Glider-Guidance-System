% Analyze MHD:
% Calculated MHD for the LC contours in mhd_osse_hindcasts_hycom_nemoV1.m
% LC is extracted in mhd_LCLE_OSSEhndcst_hycom_nemo.m
%
% NEMO LC is extracted in hycom_TSIS/extr_lc_ssh_nemo.m
%
% 
% Plot time series:
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

IEXX = zeros(9,1);  % hindcast free run # expt
IEXX(2:9) = 1;
%ixx = 3;  % hindcast/free run # expt

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

% HYCOM-TSIS hindcast experiments:
fhnd = 'hycom_tsis_expts.mat';
load(fhnd);
%Nruns = length(EXPT);
Nruns = 9; % >9 experiments are Predcitability 1a and 1b

for ii=1:Nruns
  if IEXX(ii)==0
    fprintf('%i : OFF    %s \n',ii,EXPT(ii).Name);
  else
    fprintf('%i : ON ---> %s \n',ii,EXPT(ii).Name);
  end
end

CLR = [0.4 0.8 0.2;...
       0 0.6 0.9; ...
       0 1 0.7; ...
       0 0.4 0.5; ... 
       0. 0.7 0.3; ...
       1 0 0.9; ...
       0.8 0 0.4; ...
       1 0.8 0; ...
       0.8 0.5 0; ...
       0.6 0.6 0.6];

icc=0;
ymx = 0;
for ixx=1:Nruns
  nmexp = EXPT(ixx).Name;

  if IEXX(ixx)==0, continue; end;
  fmat1 = sprintf('%sMHD_LCLCE_OSSEhndcst_hycom%2.2i_vs_nemo.mat',pthmat,ixx);
  fprintf('Loading %s\n',fmat1);
  A=load(fmat1);
  
  icc=icc+1;
  MHD(ixx).Name = EXPT(ixx).Name;
%  MHD(ixx).TM   = A.TMN;
  MHD(ixx).mhd  = A.MHD(:,1);
  mhd = MHD(ixx).mhd; 
  ymx = max([ymx, max(mhd)]);
  
end
 

TM=A.TMN;
Td=TM-TM(1)+1;
Td=Td(:);
DV=datevec(TM);

cc=0;
for iyr=2011:2012
  for im=1:12
    i1=min(find(DV(:,1)==iyr & DV(:,2)==im));
    if ~isempty(i1),
      cc=cc+1;
      ttck(cc,1)=i1;
      tlbl{cc}=sprintf('%2.2i',im);
    end
  end
end


%fnrm = 1/ymx;  % normalize scores
fnrm = 1;
ylm2 = 1.1*ymx;
if fnrm<1.0
		for ixx=1:Nruns
				if IEXX(ixx)==0; continue; end;
				MHD(ixx).mhd = MHD(ixx).mhd*fnrm; 
		end
  ylm2=1.1;
end

% ----------------------
%
%  Plot time series of MHD
%
% -----------------------

figure(1); clf;
axes('Position',[0.08 0.45 0.85 0.47]);
hold on;
for ixx=1:Nruns
  if IEXX(ixx)==0; continue; end;
  clr = CLR(ixx,:);
  mhd = MHD(ixx).mhd; 
  plot(mhd,'-','Color',clr,'Linewidth',2);
end
nrc = length(mhd);
set(gca,'tickdir','out',...
        'xlim',[1 nrc],...
        'ylim',[0 ylm2],...
        'xtick',ttck,...
        'xticklabel',tlbl,...
        'xgrid','on',...
        'ygrid','on');

if fnrm==1
  title('MHD (km) HYCOM vs NEMO LC contours, 2011/2012');
else
  title('Normalized MHD  HYCOM vs NEMO LC contours, 2011/2012');
end

%lgd = legend('HYCOM osse0','HYCOM osseE','HYCOM free','HYCOM AllSat');
%set(lgd,'position',[0.75 0.8 0.2 0.14]);
xlabel('Months');

% Legend:
axes('Position',[0.05 0.03 0.4 0.3]);
hold;
x0=0;
y0=Nruns;
icc=0;
for ixx=1:Nruns
  if IEXX(ixx)==0; continue; end;
  clr = CLR(ixx,:);
  nm = MHD(ixx).Name;
  plot([0 0.5],[y0 y0],'-','Color',clr,'Linewidth',2);
  text(0.6,y0,nm,'Fontsize',12);
  icc=icc+1;
  y0=Nruns-icc;
end
set(gca,'xlim',[-0.2 4], ...
        'ylim',[-0.1 Nruns+0.2],...
        'xtick',[],...
        'ytick',[],...
        'visible','off');

btx = 'anls_mhd_LCLCE_OSSEhndcst.m'; 

bottom_text(btx,'pwd',1,'Position',[0.52 0.04 0.4 0.04]);


%
%
ylc2=0;
SC=[];
for ixx=1:Nruns
  if IEXX(ixx)==0; continue; end;
  clr = CLR(ixx,:);
  mhd = MHD(ixx).mhd;
%
% Mean scores and StError 95%CI
  m0=nanmedian(mhd);
  Mmhd(ixx,1) = m0;
	fc1 = prctile(mhd,25);
	fc2 = prctile(mhd,75);
  Mmhd(ixx,2) = fc1;
  Mmhd(ixx,3) = fc2; 
end



% Rank by median MHD
[BB,I] = sort(Mmhd(:,1),'ascend');
axes('Position',[0.5 0.05 0.4 0.3]);
hold;
for ixx=1:Nruns
  ii0=I(ixx);
  clr=CLR(ii0,:);
  dx=0.3;
  dy=BB(ixx);
  patch([ixx-dx ixx-dx ixx+dx ixx+dx],[0 dy dy 0],clr,'Edgecolor','none');
end
set(gca,'tickdir','out',...
        'xtick',[1:10],...
        'xlim',[0 Nruns]);
ylabel('Cumulative MHD');
title('Ranking of the hindcasts');


% Legend:
axes('Position',[0.05 0.03 0.4 0.3]);
hold;
x0=0;
y0=Nruns;
icc=0;
for ixx=1:Nruns
  if IEXX(ixx)==0; continue; end;
  clr = CLR(ixx,:);
  nm = MHD(ixx).Name;
  plot([0 0.5],[y0 y0],'-','Color',clr,'Linewidth',2);
  text(0.6,y0,nm,'Fontsize',12);
  icc=icc+1;
  y0=Nruns-icc;
end
set(gca,'xlim',[-0.2 4], ...
        'ylim',[-0.1 Nruns+0.2],...
        'xtick',[],...
        'ytick',[],...
        'visible','off');

bottom_text(btx,'pwd',1,'Position',[0.002 0.04 0.4 0.04]);



% ----------------------------
% Plot Mean and St. Error: uncertainty of the mean estimate
% serr = sgm/sqrt(n)
% ----------------------------
% Rank by median MHD
Mmhd(Mmhd==0)=nan;
[MM,J] = sort(Mmhd(:,1),'ascend');


figure(4); clf;
set(gcf,'Position',[1343         734        1150         594]);
axes('Position',[0.1 0.5 0.85 0.42]);
hold;
for ixx=1:Nruns
  ii0=J(ixx);
  clr=CLR(ii0,:);
  dx=0.3;
  dy=MM(ixx,1);
  patch([ixx-dx ixx-dx ixx+dx ixx+dx],[0 dy dy 0],clr,'Edgecolor','none');

  llw = Mmhd(ii0,2);
  lup = Mmhd(ii0,3);
  plot([ixx-0.05 ixx+0.05],[llw llw],'k-');
  plot([ixx-0.05 ixx+0.05],[lup lup],'k-');
  plot([ixx ixx],[llw lup],'k-');
end
set(gca,'tickdir','out',...
        'xtick',[1:10],...
        'ygrid','on',...
        'xlim',[0 Nruns]);
if fnrm==1
  ylabel('MHD km');
else
  ylabel('Normalized MHD Score');
end
xlabel('Ranks');
title('Ranking of the hindcasts (median & IQR MHD ), LC/LCE contour');

% Legend:
axes('Position',[0.05 0.03 0.4 0.35]);
hold;
x0=0;
y0=Nruns;
icc=0;
for ixx=1:Nruns
  if IEXX(ixx)==0; continue; end;
  clr = CLR(ixx,:);
  nm = MHD(ixx).Name;
  plot([0 0.5],[y0 y0],'-','Color',clr,'Linewidth',8);
  text(0.6,y0,nm,'Fontsize',12);
  icc=icc+1;
  y0=Nruns-icc;
end
set(gca,'xlim',[-0.2 4], ...
        'ylim',[-0.1 Nruns+0.2],...
        'xtick',[],...
        'ytick',[],...
        'visible','off');

bottom_text(btx,'pwd',1,'Position',[0.002 0.04 0.4 0.04]);





