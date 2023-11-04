% Predictability experiments with
% atmospheric perturbed experiments
%
% experiments 1c - compare hycom perturbed f/casts 
% to the hycom control f/cast
% contolr run is done for 1a experiments (initialized from 
%  interpolated NEMO+GLORYS fields run for 100 days)
%
%  period =1 (start in 2011 May - June), NR - run # (01 - starts
%  on day 1 of the time period, then 7-day shift)
%
% Plot SSH and LC LCE contours used for calculating MHD
% contours are plotted for all runs within the iFcst = 10, ..., 16
% see mhd_LCLCEcntr_nemo_hycom.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;
startup;

close all
clear

% ----------------
% Flags
% ---------------

s_fig=0;

%dplot = datenum(2011,05,15); % plot date
% model day wrt to day 1 of the control run that matches the perturbed runs
IRC = [22; 52; 82];
nirc = length(IRC); 
iFcst = 11;  % 10,..., 16
% Specify perturbation forecast run: irun=2,..,5 
itime = 2;  % =1 2011, =2 2012
irun1 = 2;  % f/cast # within the time period - 7-day shifted
irun2 = 5;

dDay = [0;-2; -1; 1; 2];

% Hindcasts:
%pthd1 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
%pthd2 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
%pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
pthfrm = '/Net/kronos/ddmitry/hycom/TSIS/FIG/frames_fcst2/';

btx='plot_sshPrdct1c_LCLCEcntr.m';

% HYCOM-TSIS hindcast experiments:
fhnd = 'hycom_tsis_expts.mat';
load(fhnd);
Nruns = length(EXPT);

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

% 
% HYCOM output f/cast:
FCST = sub_fcstPrdct_info(iFcst);

irun = 1; % control run
pthhycom = FCST.TIME0(itime).RUN(irun).pthbin; % output dir
TM0 = FCST.TIME0(itime).RUN(irun).TM;

% GoM region HYCOM:
% similar to NEMO GOM domain
% for getting same endpoints of the contours
GOM=[366   489
   476   531
   542   560
   576   646
   508   827
   336   848
   204   829
    64   798
    19   746
    16   662
    12   578
    25   455
    71   382
   165   356
   281   400];

[XM,YM]=meshgrid([1:n],[1:m]);
INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
clear XM YM


cmp=flipud(colormap_cold(360));

for ii=1:nirc
  irc0 = IRC(ii);

		figure(ii); clf;
		set(gcf,'Position',[1364  550 1168 792]);

		%irc0 = find(TM0==dplot);
		dnmb0= TM0(irc0);
		DV   = datevec(dnmb0);
		yr   = DV(1);
		mo   = DV(2);
		dm   = DV(3);
		iday = dnmb0-datenum(yr,1,1)+1;

		% HYCOM   
		% Get HYCOM ssh:
		fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthhycom,yr,iday);
		finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthhycom,yr,iday);
		fprintf('Reading %s\n',fina);

		huge  = 2e20;
		rg    = 9806;

		fld = 'srfhgt';
		[F,nn,mm,ll] = read_hycom(fina,finb,fld);
		F(F>huge)=nan;
		ssh=squeeze(F)./(1e-3*rg);  % ssh m
		dmm=ssh;
		dmm(INH==0)=nan;
		%  dmm(HH>-200)=nan;
		sshM=nanmean(nanmean(dmm));
		ssh=ssh-sshM;

		%
		% Get contours from HYCOM perturbation epxeriments (iruns=2,..,5)
		nmexp = EXPT(iFcst).Name;
		pthd1 = EXPT(iFcst).path;


		%
		% PLOTTING
		clf;
		%  axes('Position',[0.08 0.55 0.8 0.38]);
		axes('Position',[0.09 0.22 0.86 0.7]);

		pcolor(LON,LAT,ssh); shading flat;

		colormap(cmp);
		caxis([-0.5 0.5]);
		hold;

		% Control run contour:
		% Get LC/LCE contours from the f/cast
		fmatout = sprintf('%shycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',pthmat2,iFcst,itime,1);
		fprintf('Loading %s\n',fmatout);
		load(fmatout);

		xh1 = LCXY.XY(irc0).X;
		yh1 = LCXY.XY(irc0).Y;
		lnW = -100;
		lce1 = 1;
		[Xhc0,Yhc0] = sub_combineLCLCE(xh1,yh1,LCE,lnW,irc0,lce1);

		CLR = [0 0 0; ...
									1  0.2 0.2; ...
									0.8 0.6 0; ...
									0.9 0 1; ...
									0.2 1 0.]; 

		for irun=irun1:irun2
		%
		% Get LC/LCE contours from the f/cast
				fmatout = sprintf('%shycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',pthmat2,iFcst,itime,irun);
				fprintf('Loading %s\n',fmatout);
				load(fmatout);

    day_shift = dDay(irun);
				TM  = LCXY.TM;
				irc  = find(TM==dnmb0-day_shift);  % match date of the control run

				xh1 = LCXY.XY(irc).X;
				yh1 = LCXY.XY(irc).Y;
				lnW = -100;
				lce1 = 1;
				[Xhc,Yhc] = sub_combineLCLCE(xh1,yh1,LCE,lnW,irc,lce1);

				plot(Xhc0,Yhc0,'.','Color',[0 0 0]);
				clr = CLR(irun,:);
				plot(Xhc,Yhc,'.','Color',clr);

		%		contour(LON,LAT,HH,[0 0],'k');

		end

		axis('equal');
		set(gca,'xlim',[-97.8 -80.7],...
			'ylim',[18.5 31],...
			'xtick',[-98:2:-82],...
			'ytick',[18:2:32],...
			'tickdir','out',...
			'Color',[0.85 0.85 0.85],...
			'Fontsize',12);

		stl=sprintf('ssh, HYCOM Fcst%2.2i %2.2i/%2.2i/%4.4i',iFcst,DV(3),DV(2),DV(1));
		title(stl);

		hb=colorbar('SouthOutside');
		set(hb,'Position',[0.2 0.1 0.5 0.015],...
		'Ticks',[-1:0.25:1],...
		'Fontsize',13);

		% Legend:
		axes('Position',[0.75 0.05 0.15 0.12]);
		hold on;
  clr = CLR(1,:); % control run
		x1=0;
		x2=0.1;
		y1=4;
		plot([x1 x2],[y1 y1],'-','Color',clr,'Linewidth',2);
		sxx = 'Control';
		text(x2+0.03,y1,sxx,'Fontsize',12);

		for irun=irun1:irun2
				clr = CLR(irun,:);
				y1=5-irun;
				plot([x1 x2],[y1 y1],'-','Color',clr,'Linewidth',2);
				sxx = sprintf('Days=%3i days',IRC(ii)-7-dDay(irun));
				text(x2+0.03,y1,sxx,'Fontsize',12);
		end
		set(gca,'xlim',[-0.02 0.42],...
										'ylim',[-0.1 4.5],...
										'xtick',[],...
										'ytick',[],...
										'visible','off');

		bottom_text(btx,'pwd',1,'Position',[0.08 0.05 0.4 0.04]);

end

				


    
    
    
    
    



