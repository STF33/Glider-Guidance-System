% Analyze MHD for the forecasts
% predictability experiments
% The forecast runs are initialized from the
% interpolated NEMO+GLORYS fields (vs data assimilative runs in OSSEs)
%
% Before analysis:
% extract LC/LCE contour
% compute MHD
% 
% Calculated MHD for the LC contours in 
% mhd_LCLCEcntr_nemo_fcsthycom.m
% LC is extracted in hycom_TSIS/extr_lc_hycom_nemoV1.m
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

% Figure to plot:
f_mhd_tser = 0; % plot MHD time series 


pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
btx = 'anls_fcstPrdct_mhd_LCLCEcntr.m';

iFcst1  = 10;  % Forecast to analyze
iFcst2  = 16;
Nruns  = 1;  % # of runs in 1 forecast series
Ntimes = 2;  % time windows when fcst is initialized: 1 - May/June 2011, 2- Jan/Feb 2012

%if iFcst1<10
%  pthmat = pthmat1;
%else
pthmat = pthmat2;
%end

ixx=0;
ymx = 0;
irun1=1;
irun2=Nruns;
for ifcst=iFcst1:iFcst2
  nhnd = ifcst;
  for itime=1:Ntimes  % 2011 and 2012 time windows
    for irun=irun1:irun2
      fmat1 = sprintf('%sMHD_LCLCE_nemo_persist_hycomfcst%2.2i-%2.2i%2.2i.mat',...
                  pthmat,nhnd,itime,irun);
      fprintf('Loading %s\n',fmat1);
      A=load(fmat1);  % MHD for f/cast and persistence
      DV=datevec(A.TM);
      ixx=ixx+1;
      MHD(ixx).mhd = A.MHD(:,1);
      MHD(ixx).Name = sprintf('fcst%2.2i-%2.2i%2.2i',nhnd,itime,irun);
      MHD(ixx).Date_str = sprintf('%2.2i/%2.2i/%4.4i',DV(1,3:-1:1));
      MHD(ixx).nemo0 = A.MHD(:,2); % NEMO true contour on day 1
      MHD(ixx).hycom0= A.MHD(:,3); % HYCOM contour on day 1
      MHD(ixx).iFcst=ifcst;    % Forecast group # (groupped by hindcast initial fields)
      MHD(ixx).Hindcast_Number=nhnd;
      MHD(ixx).Time_period=itime;
      MHD(ixx).Run_number=irun;
    end
  end
end
Ntot = ixx;



CLR = [0.2 0.9 0.2;...
       0   0.6 0.9; ...
       0.5   1   0.7; ...
       1   0.4 0.5; ... 
       0.  1  0.3; ...
       1   0   0.9; ...
       0.8 0   0.4; ...
       1   0.8 0; ...
       0.8 0.5 0; ...
       0.7 0.6 0.4];

% Persistence:
ipst = Ntot+1;
CLR(ipst,:) = [0.4 0.4 0.4];


%MXD = A.MAXD;
TM=A.TM;
Td=TM-TM(1)+1;
Td=Td(:);
DV=datevec(TM);


fnrm = 1;  % actual units - km
if fnrm<1.0
		for ixx=1:Ntot
				if IEXX(ixx)==0; continue; end;
				MHD(ixx).mhd = MHD(ixx).mhd*fnrm; 
		end
end

% ----------------------
%
%  Plot time series of MHD
%
% -----------------------
ifn = 0;
ixx = 0;
for ifcst=iFcst1:iFcst2
  for itime=1:Ntimes  % 2011 and 2012 time windows
    for irun=irun1:irun2
      ifn=ifn+1;
      itot = (ifcst-10+itime-1)*Nruns+irun;
      ixx = ixx+1;

      if f_mhd_tser==1
        fprintf('Plotting MHD time series \n');
        sub_plot_mhdPrdct_tser(MHD,CLR,irun,ifn,ipst,itot,ixx);
						  bottom_text(btx,'pwd',1,'Position',[0.09 0.02 0.4 0.04]);
      end

%
% Calculate MHD cum score
      mhd=MHD(itot).mhd;
      mhdP=MHD(itot).hycom0;
      MHD(itot).cumMHD=cumsum(mhd);
      MHD(itot).cumPST=cumsum(mhdP);
%     keyboard
    end  % runs
  end
end

% Mean MHD by hindcast/time groups
FCST = [10:16];
nFcst = length(FCST);
nll = length(MHD);
cMHD=[];
cMHDp=[];
icc=0;
for ill=1:nll
  iFcst = MHD(ill).iFcst;
  itime = MHD(ill).Time_period;
  ixx = find(FCST==iFcst);
  icc = (itime-1)*nFcst+ixx;
		% HYCOM mean MHD by 30 days
		dmm=MHD(ill).mhd;
		cMHD(icc,1)=nanmean(dmm(1:30));
		cMHD(icc,2)=nanmean(dmm(31:60));
		cMHD(icc,3)=nanmean(dmm(61:end));
		%
		% persist MHD by 30 days
		dmm=MHD(ill).hycom0;
		cMHDp(icc,1)=nanmean(dmm(1:30));
		cMHDp(icc,2)=nanmean(dmm(31:60));
		cMHDp(icc,3)=nanmean(dmm(61:end));
end




%-------------------------------------------------------------------------------- 
% Plot Bar diagrams for all f/casts by time groups


CLR = [0.5 0.2 1;
       1 0.2 0.5;
       0.7 0.7 0.7];

figure(40); clf;
set(gcf,'Position',[1316  805  858  508]);

for itime=1:Ntimes % summary for 2 time period of forecasts
  if itime==1
    axes('Position',[0.09 0.6 0.7 0.32]);
    j1=1;
    j2=nFcst;
  else
    axes('Position',[0.09 0.1 0.7 0.32]);
    j1=nFcst+1;
    j2=2*nFcst;
  end
		hold on;
		for jm=1:3
				cfc = cMHD(j1:j2,jm);
				cps = cMHDp(j1:j2,jm);

% Plot min/max for all experiments: range of mean MHD 
				mFc = mean(cfc);
				Fc1 = min(cfc);
				Fc2 = max(cfc);

				mPs = mean(cps);
				Ps1 = min(cps);
				Ps2 = max(cps);

				clr1=CLR(itime,:);
				dx=0.1;
				ixx=jm-dx;
				dy=mFc;
				patch([ixx-dx ixx-dx ixx+dx ixx+dx],[0 dy dy 0],clr1,'Edgecolor','none');
		% Range:
				llw = Fc1;
				lup = Fc2;
				plot([ixx-0.025 ixx+0.025],[llw llw],'k-');
				plot([ixx-0.025 ixx+0.025],[lup lup],'k-');
				plot([ixx ixx],[llw lup],'k-');

		% Persist
				clr2=CLR(end,:);
				ixx=jm+dx;
				dy=mPs;
				patch([ixx-dx ixx-dx ixx+dx ixx+dx],[0 dy dy 0],clr2,'Edgecolor','none');
		% Range:
				llw = Ps1;
				lup = Ps2;
				plot([ixx-0.025 ixx+0.025],[llw llw],'k-');
				plot([ixx-0.025 ixx+0.025],[lup lup],'k-');
				plot([ixx ixx],[llw lup],'k-');
		end
			
		set(gca,'tickdir','out',...
										'xlim',[0.5 3.5],...
										'xtick',[1:3],...
										'xticklabel',{'0-30','31-60','61-91'},...
										'Fontsize',14);
  dstr1=MHD(j1).Date_str;
  dstr2=MHD(j2).Date_str;
		stl=sprintf('Prdctb: MHD HYCOM fcst/prst, Ini:%s-%s',dstr1,dstr2);
  title(stl);
end
axes('Position',[0.8 0.6 0.15 0.25]);
hold on
clr1 = CLR(1,:);
clr2 = CLR(2,:);
clr3 = CLR(3,:);
patch([1 1. 1.2 1.2],[1 1.2 1.2 1],clr1,'Edgecolor','none');
patch([1 1. 1.2 1.2],[0.7 0.9 0.9 0.7],clr2,'Edgecolor','none');
patch([1 1. 1.2 1.2],[0.4 0.6 0.6 0.4],clr3,'Edgecolor','none');

text(1.3,1.1,'Fcst 2011','Fontsize',12);
text(1.3,0.8,'Fcst 2012','Fontsize',12);
text(1.3,0.5,'Persist','Fontsize',12);
set(gca,'xlim',[1 2],...
        'ylim',[0.3 1.2],...
        'visible','off');

bottom_text(btx,'pwd',1);		

	
% -------------------- 
%  Plot Time 1 vs Time 2 
% ----------------------------------------  
dx0=0.1;

figure(41); clf;
set(gcf,'Position',[1316  805  858  508]);
axes('Position',[0.09 0.6 0.7 0.32]);
for itime=1:Ntimes % summary for 2 time period of forecasts
  if itime==1
    j1=1;
    j2=nFcst;
    dx=-dx0;
  else
    j1=nFcst+1;
    j2=2*nFcst;
    dx=dx0;
  end
  hold on;
  for jm=1:3
    cfc = cMHD(j1:j2,jm);
    cps = cMHDp(j1:j2,jm);

% Plot min/max for all experiments: range of mean MHD 
    mFc = mean(cfc);
    Fc1 = min(cfc);
    Fc2 = max(cfc);


    clr1=CLR(itime,:);
    ixx=jm+dx;
    dy=mFc;
    patch([ixx-dx ixx-dx ixx+dx ixx+dx],[0 dy dy 0],clr1,'Edgecolor','none');
  % Range:
    llw = Fc1;
    lup = Fc2;
    plot([ixx-0.025 ixx+0.025],[llw llw],'k-');
    plot([ixx-0.025 ixx+0.025],[lup lup],'k-');
    plot([ixx ixx],[llw lup],'k-');

  end

  set(gca,'tickdir','out',...
          'xlim',[0.5 3.5],...
          'xtick',[1:3],...
          'xticklabel',{'0-30','31-60','61-91'},...
          'Fontsize',14);
  ds1=MHD(j1).Name;
  ds2=MHD(j2).Name;
  dstr1=MHD(j1).Date_str;
  dstr2=MHD(j2).Date_str;
  stl=sprintf('Predictb Expts, MHD HYCOM fcst 2011/2012');
  title(stl);
end

clr1 = CLR(1,:);
clr2 = CLR(2,:);
axes('Position',[0.8 0.7 0.15 0.15]);
hold on
patch([1 1. 1.2 1.2],[1 1.2 1.2 1],clr1,'Edgecolor','none');
patch([1 1. 1.2 1.2],[0.7 0.9 0.9 0.7],clr2,'Edgecolor','none');
text(1.3,1.1,'Fcst 2011','Fontsize',12);
text(1.3,0.8,'Fcst 2012','Fontsize',12);
set(gca,'xlim',[1 2],...
        'ylim',[0.6 1.2],...
        'visible','off');

bottom_text(btx,'pwd',1,'Position',[0.03 0.4 0.5 0.05]);






