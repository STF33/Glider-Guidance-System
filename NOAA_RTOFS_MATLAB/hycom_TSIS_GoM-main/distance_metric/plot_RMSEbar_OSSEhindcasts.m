% Plot RMSE for ssh from OSSE hindcasts
% that use NEMO as NR to assimilat synthetic observations
% Bar diagrams = median from the overall RMSE for each hindcast
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

% Set flags for extracting experiments:
EXON = zeros(9,1);
EXON(2:9) = 1; % select expt to be extracted,  #2 - ssh ???
f_nemo = 0;    % obsolete - use extr_lc_ssh_nemo.m
Bisol = 0.17;  % ssh contour

% Hindcasts:
pthd1  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
pthd2  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
pthd12 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_obs/';  % all satellites included
pthd3  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/2011_GLfreerun/'; % free run
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

%fmatout = sprintf('%sLC_coord_osse_hycomV1.mat',pthmat);

btx = 'plot_RMSEbar_OSSEhindcasts.m';

fhnd = 'hycom_tsis_expts.mat';
fprintf('Loading %s\n',fhnd);
load(fhnd);

Nruns = length(EXON);

for ii=1:Nruns
  if EXON(ii)==0;
    fprintf('%i : OFF    %s \n',ii,EXPT(ii).Name);
  else
    fprintf('%i : ON ---> %s \n',ii,EXPT(ii).Name);
  end
end


for ii=1:Nruns
 fprintf('%i: %s \n',ii,EXPT(ii).path);
end

%Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
HH(isnan(HH))=100;

Iocn = find(HH<-200); % deep GoM

SC=[];
POOL = struct;
Iexpt = find(EXON==1);
for jj=1:length(Iexpt);
  ixx = Iexpt(jj);
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;
  fmatout = sprintf('%shycom_nemo_RMSE_hcst%2.2i.mat',pthmat,ixx);
  fprintf('%s Loading %s\n',nmexp,fmatout);
  load(fmatout);

  nmexp= RMSE.Name;
  rmse = RMSE.rmse;
  TM   = RMSE.TM;
  dve=datevec(TM(end));

  dmm = rmse(Iocn);

  POOL(jj).Name=nmexp;
  POOL(jj).rmse_md = nanmedian(dmm);
  POOL(jj).rmse_p25 = prctile(dmm,25);
  POOL(jj).rmse_p75 = prctile(dmm,75);
  POOL(jj).nmb_expt = ixx;
  SC(jj) = nanmedian(dmm);
end

% Hindcast colors 
% match with hindcast for
% comparison
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

% Rank by cum score:
SC(SC==0)=nan;
[BB,I] = sort(SC,'ascend');

figure(1); clf;
set(gcf,'Position',[1343         734        1150         594]);
axes('Position',[0.1 0.4 0.8 0.4]);
hold;
for ixx=1:length(I)
  jj=I(ixx);
  enmb = POOL(jj).nmb_expt;
  clr = CLR(enmb,:);
  dx=0.3;
  dy=BB(ixx);
  patch([ixx-dx ixx-dx ixx+dx ixx+dx],[0 dy dy 0],clr,'Edgecolor','none');
%
% IQR:
  llw = POOL(jj).rmse_p25;
  lup = POOL(jj).rmse_p75;
  plot([ixx-0.05 ixx+0.05],[llw llw],'k-','linewidth',1.6);
  plot([ixx-0.05 ixx+0.05],[lup lup],'k-','linewidth',1.6);
  plot([ixx ixx],[llw lup],'k-','linewidth',1.6);

end
set(gca,'tickdir','out',...
        'xtick',[1:10],...
        'ygrid','on',...
        'xlim',[0 Nruns]);
  
ylabel('RMSE SSH (m)');

title('Median & IQR RMSE(SSH) OSSE hndcst vs NEMO');

% Legend:
axes('Position',[0.05 0.03 0.4 0.3]);
hold;
x0=0;
y0=Nruns;
icc=0;
for ixx=1:Nruns
  if EXON(ixx)==0; continue; end;
  icc=icc+1;
  clr = CLR(ixx,:);
  nm = POOL(icc).Name;
  plot([0 0.5],[y0 y0],'-','Color',clr,'Linewidth',2);
  text(0.6,y0,nm,'Fontsize',12);
  y0=Nruns-icc;
end
set(gca,'xlim',[-0.2 4], ...
        'ylim',[-0.1 Nruns+0.2],...
        'xtick',[],...
        'ytick',[],...
        'visible','off');

bottom_text(btx,'pwd',1,'Position',[0.002 0.04 0.4 0.04]);




