% Calculate MHD for the LC 
% in the Eastern GoM
% East of 90W
%
% for the HYCOM forecasts and NEMO
% for reference - persistence: day 0 of the free run
% is used for 2011 only
% extracted in hycom_TSIS/distance_metric/extr_lc_hycom_fcst.m
%
%  Forecasts: decision tree:
%
%  Hindcast group (which is used to initialize f/cst): <10 - OSSEs, 10 - predictab f/casts
%      |
%      V
%  TimePeriod during which the f/casts are initialized (May-June 2011, Jan-Feb 2012)
%     |
%     V
%  Forecasts (90 day f/csts shifted 7 days - 7 fcsts)
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear


% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat0  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';

Nruns = 1; % # of forecast runs - hindcast groups
NTime = 1; % Time periods of the forecasts: started in May/June 2011 and Jan/Feb 2012 
Nfcst = 7; 

% All hindcast experiments
load('hycom_tsis_expts.mat');

% Array of forecast runs
icc=0;
irun1=1;
irun2=Nfcst;
for ifcst=10:10
 switch(ifcst)
  case(1)
    Nhnd=7;
  case(2)
    Nhnd=8;
  case(3)
    Nhnd=2;
    NTime=1;
    irun1=2;
    irun2=2;
   case(4)
    Nhnd=3;
   case(10);
    Nhnd=10;
    NTime=1;
    irun1=1;
    irun2=1;
  end
  for it=1:NTime
    for irun=irun1:irun2
      icc=icc+1;
      FRUNS(icc).Nhndcst = Nhnd;
      FRUNS(icc).TimePeriod = it;
      FRUNS(icc).Forecast_group = ifcst;
      FRUNS(icc).Nfcst = irun;
      FRUNS(icc).matname = sprintf('hycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',Nhnd,it,irun);
      FRUNS(icc).Hind_name = EXPT(Nhnd).Name;
    end
  end
end

fmat = 'hycom_info_forecasts.mat';
fprintf('saving forecast info %s\n',fmat);
save(fmat,'FRUNS');


Mruns = icc;

% YC point:
x0  = -85;
y0  =  22;
lnW = -90; % cutoff longitude

%if f_mhd==1
% Array with NEMO LC
%  fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat); OLD
fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat1);
fprintf('Loading %s\n',fmatout);
load(fmatout);
LCN=LCXY;
LCEN=LCE;  % LCE contours
TMN = LCN(1).TM;  % Nemo

clear LCXY LCE

for ixx=1:Mruns
		nmexp = FRUNS(ixx).Hind_name;

  flnm = FRUNS(ixx).matname;
  Nhndcst = FRUNS(ixx).Nhndcst; 
  if Nhndcst>=10
    pthmat = pthmat2;
  else
    pthmat = pthmat0;
  end
		fmatout = sprintf('%s%s',pthmat,flnm);
		fprintf('Loading %s\n',fmatout);
		load(fmatout);

		TM  = LCXY.TM;
  nrc = length(TM);
% Persistence contours: initial state NEMO
  xPs = [];
  yPs = [];
  xhPs= [];
  yhPs= [];
  MHD = [];

		fprintf('Calculating MHD %s\n',nmexp);
		for ihc=1:nrc
				if mod(ihc,30)==0,
						fprintf('  %6.2f%% done ...\n',ihc/nrc*100);
				end

				dm1= TM(ihc);
				irc= find(TMN==dm1);

				Xn = LCN(1).XY(irc).X; % NEMO LC contour
				Yn = LCN(1).XY(irc).Y;

%    [Xc,Yc] = sub_combineLCLCE(xn,yn,LCEN,lnW,irc);

% Save 1st contour NEMO - persistence reference  
    if ihc==1
      xPs=Xn;
      yPs=Yn;
    end
   
				if ~isempty(irc);
						Xhc = LCXY.XY(ihc).X;
						Yhc = LCXY.XY(ihc).Y;

%     [Xhc,Yhc] = sub_combineLCLCE(xh1,yh1,LCE,lnW,ihc);

% Save HYCOM 1st contour for persistence
      if ihc==1
        xhPs = Xhc;
        yhPs = Yhc;
      end

						P = [Xn,Yn];    % NEMO
						Q = [Xhc,Yhc];  % HYCOM fcst
						mhd1 = modified_hausdorff_distance(P,Q,'geo');

%      if mhd1>100
%        fprintf('mhd1>100 ixx=%i, ihc=%i %s\n',ixx,ihc,datestr(dm1));
%        keyboard
%      end

% Persistence using NEMO
      Q = [xPs,yPs];
      mhdPS = modified_hausdorff_distance(P,Q,'geo');

% Persistence - 1st day of HYCOM - not identically same as NEMO
      Q = [xhPs,yhPs];
      mhdPS_hycom = modified_hausdorff_distance(P,Q,'geo');

				else
						mhd1 = nan;
      mhdPS = nan;
      mhdPS_hycom = nan;
				end;

				MHD(ihc,1)  = mhd1;         % mhd fcst - nemo
    MHD(ihc,2)  = mhdPS;        % persistence nemo - nemo
    MHD(ihc,3)  = mhdPS_hycom;  % persistence hycom - nemo
    
%keyboard
		end
%keyboard				
		fprintf('LC: min/max MHD=%8.2f  %8.2f, min/max prst MHD=%8.2f %8.2f\n',...
          min(MHD(:,1)),max(MHD(:,1)),min(MHD(:,2)),max(MHD(:,2)));

		fmat1 = sprintf('%sMHD_LC_nemo_persist_hycomfcst%2.2i-%2.2i%2.2i.mat',...
                  pthmat,FRUNS(ixx).Nhndcst,FRUNS(ixx).TimePeriod,FRUNS(ixx).Nfcst);
		fprintf('Saving %s\n',fmat1);
		save(fmat1,'MHD','TM');
%keyboard	
end;
%end

fprintf('All done\n');

f_plt=0;
if f_plt==1
  figure(1); clf;
  axes('Position',[0.09 0.5 0.84 0.42]);
  hold on;
  plot(MHD(:,1),'linewidth',2);
  plot(MHD(:,2),'linewidth',2);

  lgd = legend('HYCOM-NEMO','Persist');
  set(lgd,'Position',[0.2 0.3 0.3 0.12]);

  set(gca,'tickdir','out',...
      'xlim',[0 90],...
      'ylim',[0 max(max(MHD(:,1:2)))],...
      'xtick',[0:10:100],...
      'xgrid','on',...
      'ygrid','on',...
      'fontsize',13);

  fnm = FRUNS(ixx).matname;
  dv1 = datevec(TM(1));
  dv2 = datevec(TM(end));
  stl = sprintf('MHD, %s, %s-%s',fnm,datestr(dv1),datestr(dv2));
  title(stl,'Interpreter','none');

  btx = 'mhd_LCcntr_nemo_fcsthycom.m';
  bottom_text(btx,'pwd',1);

end





