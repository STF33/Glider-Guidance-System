% Calculate MHD for the LC and LCE contours 
% in the Eastern GoM
% East of 90W
%
% for the HYCOM hindcasts and NEMO
% for reference - free run HYCOM
% is used for 2011 only
% extracted in hycom_TSIS/extr_lc_hycom_nemo.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

%f_mhd = 1;  % = 1 - calculate MHD, =
IEXX = zeros(9,1);
IEXX(2:9) = 1;  % hindcast free run # expt
%ixx = 3;  % hindcast/free run # expt

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

% HYCOM-TSIS hindcast experiments:
fhnd = 'hycom_tsis_expts.mat';
load(fhnd);
Nruns = length(EXPT);

for ii=1:Nruns
  if IEXX(ii)==0
    fprintf('%i : OFF    %s \n',ii,EXPT(ii).Name);
  else
    fprintf('%i : ON ---> %s \n',ii,EXPT(ii).Name);
  end
end


% YC point:
x0  = -85;
y0  =  22;
lnW = -90; % cutoff longitude

%if f_mhd==1
% Array with NEMO LC
%  fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat); OLD
fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat);
fprintf('Loading %s\n',fmatout);
load(fmatout);
LCN=LCXY;
LCEN=LCE;  % LCE contours

TMN = LCN(1).TM;  % Nemo
nrc = length(LCN(1).XY);

for ixx=1:Nruns
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;

  if IEXX(ixx)==0, continue; end; 

  fmatout = sprintf('%shycom_LCcontour_%2.2i.mat',pthmat,ixx);
  fprintf('Loading %s\n',fmatout);
  load(fmatout);

  TM  = LCXY.TM;

  fprintf('Calculating MHD %s\n',nmexp);
  for irc=1:nrc
    if mod(irc,100)==0,
      fprintf('  %6.2f%% done ...\n',irc/nrc*100);
    end

    dm1= TMN(irc);
    ihc= find(TM==dm1);

    xn = LCN(1).XY(irc).X; % NEMO LC contour
    yn = LCN(1).XY(irc).Y;

   [Xc,Yc] = sub_combineLCLCE(xn,yn,LCEN,lnW,irc);

    if ~isempty(ihc);
      xh1 = LCXY.XY(ihc).X;
      yh1 = LCXY.XY(ihc).Y;

     [Xhc,Yhc] = sub_combineLCLCE(xh1,yh1,LCE,lnW,irc);


      P = [Xc,Yc];
      Q = [Xhc,Yhc];
      mhd1 = modified_hausdorff_distance(P,Q,'geo');

    else
      mhd1 = nan;
    end;

    MHD(irc,1)  = mhd1;
%keyboard
  end
%keyboard    
  fprintf('  min/max MHD=%8.2f  %8.2f\n',min(MHD),max(MHD));
  fmat1 = sprintf('%sMHD_LCLCE_hycom%2.2i_vs_nemo.mat',pthmat,ixx);
  fprintf('Saving %s\n',fmat1);
  save(fmat1,'MHD','TMN');
 
end;
%end

fprintf('All done\n');


