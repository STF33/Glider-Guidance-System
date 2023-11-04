% Plot RMSE for ssh
% sqrt (sum[(a-<a>)^2]/n) 
%
% Extract and save LC contours from HYCOM-TSIS and
% nemo simulations with new hindcasts and free run 
% specify individually which run need to extract
%
%  NEMO is extracted in extr_lc_hycom_nemo-V0.m 
%
% Compare LC contours from the HYCOM_TSIS hindcast
% assimilating NEMO 1/100 free running fields
% and NEMO LC contours
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

if f_nemo>0
  fprintf('NEMO is now extracted in extr_lc_ssh_nemo.m !!! flag ignored\n\n');
  f_nemo=0;
end

% Hindcasts:
pthd1  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
pthd2  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
pthd12 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_obs/';  % all satellites included
pthd3  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/2011_GLfreerun/'; % free run
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

%fmatout = sprintf('%sLC_coord_osse_hycomV1.mat',pthmat);

btx = 'plot_rmse_ssh.m';


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


CMP = create_colormapBGY(400,0,0.2);
cmp = CMP.colormap;
c1=0;
c2=0.2;

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

  figure('Position',[1238 482 1280 857]);
  clf; hold on;
  pcolor(LON,LAT,rmse); shading flat;
  colormap(cmp);
  caxis([c1 c2]);
  chb = colorbar;
  

  contour(LON,LAT,HH,[0 0],'-','Color',[1 1 1]);
  contour(LON,LAT,HH,[-500 -500],'-','Color',[0.8 0.8 0.8]);


  axis('equal');
  set(gca,'tickdir','out',...
          'xlim',[-98 -79],...
          'ylim',[18 31],...
          'Color',[0 0 0], ...
          'Fontsize',14);

  stl=sprintf('RMSE SSH (hycom-nemo), %2.2ihcst: %s, %s-%s',ixx,nmexp,datestr(TM(1)),datestr(TM(end)));
  title(stl);

  bottom_text(btx,'pwd',1);




end


