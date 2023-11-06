% Calculate MHD using ALL LCEs + LC
%
% HYCOM Predictability experiments: NEMO+GLORYS interpolated for
% Experiments 1b: 
%  
% For each f/cast: 10, ..., 16
% there are: control run = iRun =1
% and  4 perturbed simulations, iRun=2,..., 5
%
% for the HYCOM forecasts and HYCOM control run
% for reference - persistence: day 0 of the free run
% is used for 2011 only
% extracted in hycom_TSIS/distance_metric/extr_lc_hycom_fcst.m
%
%
%  Predictability experiments: interpolated NEMO+GLORYS --> HYCOM
%        #10, ..., #16

%  Forecast group (which is used to initialize f/cst) iFcst
%      |
%      V
%  TimePeriod during which the f/casts are initialized (May-June 2011, Jan-Feb 2012)
%     |
%     V
%  irun: free running forecasts (90 day f/csts) = 2, ..., 5 - perturb f/casts
%        = 1 - control run (initialized from NEMO) 
%
% FSU COAPS, D. Dukhovskoy
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

% 1 forecast at a time
iFcst=14; % Prdct1b f/casts 10, ..., 16 
itime1=1;
itime2=2;
irun1=2;  % irun=2,...,5
irun2=5;

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';

btx = 'mhd_allLCLCEcntrPrdct1b_hycom0_fcsthycom.m';

% All hindcast experiments
%load('hycom_tsis_expts.mat');

FCST = sub_fcstPrdct_info(iFcst);


% YC point:
x0  = -85;
y0  =  22;

clear LCXY LCE


iTime = 0;
for itime=itime1:itime2
% Control run - HYCOM without perturbation   
  iTime = itime; 
	irun0 = 1; % control run => irun=1
	flnm = sprintf('hycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',iFcst,itime,irun0);    
	fmatcntr = sprintf('%s%s',pthmat2,flnm);
	fprintf('Loading Control run %s\n',fmatcntr);
	load(fmatcntr);
	TMN = LCXY.TM;
	LCN = LCXY;
	LCEN = LCE;

  for irun=irun1:irun2
    RUN0  = FCST.TIME0(itime).RUN(1); % control run
    RUN   = FCST.TIME0(itime).RUN(irun);
    pthd1 = RUN.pthbin;
    TM    = RUN.TM;
    YDAY  = RUN.jday;
    nrc   = length(TM);
    day_shift = RUN.dayT0_true-RUN.prdct_time0;
    TM    = TM+day_shift;
    DV    = datevec(TM);

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
    fmat1 = sprintf('%sMHD_allLCLCE_Prdct1b_%s.mat',pthmat2,fcstname);
    nmexp = fcstname;
  

    fprintf('\n\n %s %s\n',fcstname,fmat1);
    fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',iFcst,itime,irun);
    fprintf(' Input data: %s\n',pthd1);
    fprintf(' Output saved: %s\n',fmat1);

    flnm = sprintf('hycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',iFcst,iTime,irun);
    fmatfcst = sprintf('%s%s',pthmat2,flnm);
    fprintf('Loading f/cast run %s\n',fmatfcst);
    load(fmatfcst);

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

			xn = LCN(1).XY(irc).X; % control LC contour
			yn = LCN(1).XY(irc).Y;
%
% Compare f/cast with all LCEs to the control run
      Nlcec = length(LCEN);
      Nlcef = length(LCE);
%
% Calculate MHD using ALL LCEs 
      if ~isempty(irc);
				xh1 = LCXY.XY(ihc).X;
				yh1 = LCXY.XY(ihc).Y;

        MHD_LCLCE = [];
% 
% All LCEs to the LC
        Xc = xn;
        Yc = yn;
        for ilce=1:Nlcec+1   % adding eddies to the control run contour
          Neddy = ilce-1;
          [Xc,Yc,flgC] = sub_combineLCLCE_EddyN(Xc,Yc,LCEN,irc,Neddy);
        end

        Xhc=xh1;
        Yhc=yh1;
				for jlce=1:Nlcef+1   % adding eddies to the f/cast
					Neddy = jlce-1;
					[Xhc,Yhc,flg] = sub_combineLCLCE_EddyN(Xhc,Yhc,LCE,ihc,Neddy);
        end

				P = [Xc,Yc];    % HYCOMcontrol run
				Q = [Xhc,Yhc];  % HYCOM fcst LC 
				mhd_fcst = modified_hausdorff_distance(P,Q,'geo');
%
% --------------------------------------------------------
      f_chck_ssh = 89;
      if f_chck_ssh==ihc
        fprintf('Plotting SSH and contours for checking\n');
				DV1 = datevec(dm1);
        stl=sprintf('ssh, MHD=%7.4g, HYCOM fcst%i run=%2.2i, day=%3.3i,  %2.2i/%2.2i/%4.4i',...
                   mhd_fcst,iFcst,irun,ihc,DV1(3),DV1(2),DV1(1));

        sub_plot_ssh_LCLCE_chck(iFcst,itime,irun,dm1,Q,P,FCST,stl)

        bottom_text(btx,'pwd')
        keyboard
      end
% --------------------------------------------------------
%
% Same for persistence 
% Persistence using day=t0 from HYCOM control run
        if ihc==1
          xPs = Xc;
          yPs = Yc;
        end

        Q  = [xPs,yPs];
        P = [Xc,Yc];    % HYCOMcontrol run
        mhdPs = modified_hausdorff_distance(P,Q,'geo');

				if mhd_fcst>100
				  fprintf('WARNING: high  mhd %4.1f %s\n',mhd_fcst,datestr(dm1));
				end


			else
					mhd_fcst = nan;
					mhdPs = nan;
			end;

			MHD(ihc,1)  = mhd_fcst;         % mhd fcst - hycom control run
			MHD(ihc,2)  = mhdPs;            % persistence hycom (t=t0)
						
		end

		fprintf('  min/max MHD=%8.2f  %8.2f, min/max prst MHD=%8.2f %8.2f\n',...
				min(MHD(:,1)),max(MHD(:,1)),min(MHD(:,2)),max(MHD(:,2)));
%keyboard
		fprintf('Saving %s\n',fmat1);
		save(fmat1,'MHD','TM');

  end
end;
%end

fprintf('All done\n');


