% Calculate MHD for the LC and LCE contours 
% in the Eastern GoM
% East of 90W
%
% for the HYCOM hindcasts and NEMO
% for reference - free run HYCOM
% is used for 2011 only
% extracted in hycom_TSIS/extr_lc_temp 
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

T0 = 2.5;  % dt contour
Z0 = -200;

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
fTNM = sprintf('%sNEMO_dT%2.2i_Z%3.3i_contour.mat',pthmat,round(T0*10),abs(Z0));
fprintf('Loading %s\n',fTNM);
load(fTNM);
LCN=LCXY;


TMN = LCN(1).TM;  % Nemo
nrc = length(LCN(1).XY);

for ixx=1:Nruns
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;

  if IEXX(ixx)==0, continue; end; 

  fmatout2 = sprintf('%sHYCOM_dT%2.2i_Z%3.3i_contour_hind%2.2i.mat',...
                    pthmat,round(T0*10),abs(Z0),ixx);
  fprintf('Loading %s\n',fmatout2);
  load(fmatout2);

  LCH=LCXY;
  TM = LCH.TM;



  fprintf('Calculating MHD %s\n',nmexp);
  for irc=1:nrc
    if mod(irc,100)==0,
      fprintf('  %6.2f%% done ...\n',irc/nrc*100);
    end

    dm1= TMN(irc);
    ihc= find(TM==dm1);

    Xn = LCN(1).XY(irc).X; % NEMO LC contour
    Yn = LCN(1).XY(irc).Y;

% Clean contours near Cuba:
%    lnW = -90;
    xc1=-86.8;
    xc2=-84.07;
    xc3=xc2;
    xc4=-81.1;
    yc1=22.8;
    yc2=yc1;
    yc3=25.1;
    yc4=yc3;
    I=find(Yn<=yc1 & Xn>xc1);
    Xn(I)=nan;
    Yn(I)=nan;
    I=find(Yn<=yc3 & Xn>xc3);
    Xn(I)=nan;
    Yn(I)=nan;
    I=find(Xn<lnW);
    Xn(I)=nan;
    Yn(I)=nan;
    Inn = find(~isnan(Xn));
    Xn=Xn(Inn);
    Yn=Yn(Inn);

    if ~isempty(ihc);
      Xh = LCH.XY(ihc).X;
      Yh = LCH.XY(ihc).Y;

      if isempty(Xh)
        fprintf('HYCOM contour is empty \n');
        keyboard;
      end

      I=find(Yh<=yc1 & Xh>xc1);
      Xh(I)=nan;
      Yh(I)=nan;
      I=find(Yh<=yc3 & Xh>xc3);
      Xh(I)=nan;
      Yh(I)=nan;
      I=find(Xh<lnW);
      Xh(I)=nan;
      Yh(I)=nan;
      Inn = find(~isnan(Xh));
      Xh=Xh(Inn);
      Yh=Yh(Inn);



      P = [Xn,Yn];
      Q = [Xh,Yh];
      mhd1 = modified_hausdorff_distance(P,Q,'geo');

    else
      mhd1 = nan;
    end;

    if isempty(mhd1)   % no contour
      MHD(irc,1) = nan;
    else
      MHD(irc,1)  = mhd1;
    end

%keyboard
  end
%keyboard    
  fprintf('  min/max MHD=%8.2f  %8.2f\n',min(MHD),max(MHD));
  fmat1 = sprintf('%sMHD_dT%2.2i_Z%3.3i_hycom%2.2i_vs_nemo.mat',...
               pthmat,round(T0*10),abs(Z0),ixx);
  fprintf('Saving %s\n',fmat1);
  save(fmat1,'MHD','TMN');
 
end;
%end

fprintf('All done\n');


