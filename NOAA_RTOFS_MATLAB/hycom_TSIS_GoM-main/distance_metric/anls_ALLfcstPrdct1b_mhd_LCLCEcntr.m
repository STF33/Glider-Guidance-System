% Experiments 1b: HYCOM-HYCOM 
%
% PRedictability experiments:
% initial fields from HYCOM control run
% Control run - f/cast with initical fields interpolated from HYCOM+GLORYS
%
% Analyze MHD for all forecasts
% 
% Before analysis:
% extract LC/LCE contour
%  extr_lc_hycomPrdct_fcst.m
% compute MHD
% mhd_LCLCEcntrPrdct1b_hycom0_fcsthycom.m 
%
% MHD is calculated in other codes:
% There several approaches used for defining LC:
% LC + 1 LCE that gives the best MHD: mhd_LCLCEcntrPrdct1b_hycom0_fcsthycom.m 
%
% LC + all LCEs: mhd_allLCLCEcntrPrdct1b_hycom0_fcsthycom.m
%
% LC is extracted in hycom_TSIS/extr_lc_hycom_nemo.m
% NEMO LC is extracted in hycom_TSIS/extr_lc_ssh_nemo.m
%
% Control runs:
% Forecasts: #10 - May 01, 2011 - 100 days
%            #11 - May 08, 2011 
%            #12 - May 15, 2011
% ...
%            #16 - May 
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

ifcst1=10;
ifcst2=16;
irun1 = 2; % irun=2,...,5 
irun2 = 5;
itime1= 1;
itime2= 2;

% Figure to plot:
f_mhd_tser = 0; % plot daily MHD time series 

DT = [0; -2; -1; 1; 2];

%pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
btx = 'anls_ALLfcstPrdct1b_mhd_LCLCEcntr.m';


% Hindcast info
load('hycom_tsis_expts.mat');  % EXPT array

ixx=0;
ymx = 0;
ii=0;
for iFcst=ifcst1:ifcst2
  ii=ii+1;
  for itime=itime1:itime2  % 2011 and 2012 time windows
    irr=0;
    ymax = -10;
    for irun=irun1:irun2
      fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
%      fmat1 = sprintf('%sMHD_LCLCE_Prdct1b_%s.mat',pthmat,fcstname);  % best MHD LC+1LCE
      fmat1 = sprintf('%sMHD_allLCLCE_Prdct1b_%s.mat',pthmat,fcstname);   % LC + all LCEs

      fprintf('Loading %s\n',fmat1);
      A=load(fmat1);  % MHD for f/cast and persistence
      TM=A.TM;
      DV=datevec(TM);
      smm = max(A.MHD(:,1:2));
      ymax= max([smm,ymax]);

      MHD(ii).TIME(itime).RUN(irun).mhd      = A.MHD(:,1);
      MHD(ii).TIME(itime).RUN(irun).Name     = fcstname;
      MHD(ii).TIME(itime).RUN(irun).Date_str = sprintf('%2.2i/%2.2i/%4.4i',DV(1,3:-1:1));
      MHD(ii).TIME(itime).RUN(irun).mhdPrst  = A.MHD(:,2); % persistence mhd
      MHD(ii).TIME(itime).RUN(irun).iFcst    = iFcst;    % Forecast group # (groupped by hindcast initial fields)
    end
    MHD(ii).TIME(itime).Ymax=ymax;
  end
end
Nfgr = ii;

% Colors for different perturbed runs within same f/cast
CLR = [0   0    0;...
       0.8 0.2  0;...
       1   0.6  0; ...
       0   1    0.5; ... 
       0   0.5  1];

% Persistence:
clrP = [0.4 0.4 0.4];


%MXD = A.MAXD;
TM=A.TM;
Td=TM-TM(1)+1;
Td=Td(:);
DV=datevec(TM);

% ----------------------
%
%  Plot time series of MHD
%
% -----------------------
ifn = 0;
clear ixx
nRuns = length(MHD);
if f_mhd_tser==1
  fprintf(' ===== Plotting MHD Time Series =====\n\n');
  ii=0;
	for iFcst=ifcst1:ifcst2
		ii=ii+1;
		for itime=itime1:itime2  % 2011 and 2012 time windows
			ifn=ifn+1;
      for irun=irun1:irun2
        dstr = MHD(ii).TIME(itime).RUN(irun).Date_str;
        fcstname = MHD(ii).TIME(itime).RUN(irun).Name;
        clr  = CLR(irun,:);
        mhd  = MHD(ii).TIME(itime).RUN(irun).mhd;
        mhdP = MHD(ii).TIME(itime).RUN(irun).mhdPrst;
        ymax = MHD(ii).TIME(itime).Ymax;
        if irun==irun1
          stl = sprintf('MHD(km) HYCOM Frcst1b %s LCLCE cntrs, %s',fcstname,dstr);
          lgnd = -1;
						    sub_plot_mhd_tser_v2(mhd,mhdP,clr,clrP,ifn,stl,ymax,lgnd);
        else
          plot(mhd,'-','Color',clr,'Linewidth',2);
        end
      end
% Legend
			axes('Position',[0.1 0.2 0.3 0.3]);
			hold on;
			x1=0.1;
			x2=x1+0.2;
			y1=1;
			for irun=irun1:irun2
				y1=y1-0.1;
				clr=CLR(irun,:);
				plot([x1 x2],[y1 y1],'-','Color',clr,'Linewidth',2); 
				dt = DT(irun); 
				stxt = sprintf('run=%i dt=%i',irun,dt);
				text(x2+0.1,y1,stxt,'Fontsize',12);
			end
			y1=y1-0.1;
			clr=clrP;
			plot([x1 x2],[y1 y1],'-','Color',clr,'Linewidth',2);
			stxt = 'Persist'; 
			text(x2+0.1,y1,stxt,'Fontsize',12);
			set(gca,'xlim',[x1 x1+0.6],...
											'ylim',[y1-0.1 1], ...
											'Visible','off');   

			bottom_text(btx,'pwd',1,'Position',[0.09 0.02 0.4 0.04]);
		end
	end

  fprintf('Time Series plotted ...\n');
  keyboard
end


% Hindcast colors 
% match with hindcast for
% comparison
CLRH = [0.4 0.8 0.2;...
       0   0.6 0.9; ...
       0.5   1   0.7; ...
       1   0.4 0.5; ... 
       0.  1  0.3; ...
       1   0   0.9; ...
       0.8 0   0.4; ...
       1   0.8 0; ...
       0.8 0.5 0; ...
       0.7 0.6 0.4; ...
       0.5 0.2 1];


% Group by time shift -2 days, -1 day, +1 day, +2 days
% show individual f/casts 
% Weekly pools
WK=[1:7:92]; % weekly pools
nwk=length(WK)-1;

% Mean MHD by forecast-hindcast/time groups
% ----------------------------
% Plot: all F/csts in 1 group for 3 30-day time periods
% coordinates of the Bars relative to X=# of the 30-day period, 1,2,3
nruns=irun2-irun1+1; 

nfg=40;
anls_nm='MHD(km) HYCOM-Prdct1b LCLCE months';
sub_plotMHD_bars_month_Prdct1b(nfg,CLRH,MHD,Nfgr,ifcst1,ifcst2,irun1,irun2,anls_nm);
bottom_text(btx,'pwd',1);		

% Plot weekly bars
nfg=41;
anls_nm='MHD(km) HYCOM-Prdct1b LCLCE weeks';
sub_plotMHD_bars_week_Prdct1b(nfg,CLRH,MHD,Nfgr,ifcst1,ifcst2,...
                               irun1,irun2,WK,anls_nm);
bottom_text(btx,'pwd',1);








