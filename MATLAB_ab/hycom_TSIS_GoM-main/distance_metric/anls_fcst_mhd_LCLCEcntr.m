% Analyze MHD for the forecasts
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
% Hindcasts used for the forecast runs:
% H/cast #2 - Full 2D SSH                         Fcst 3
%        #3 - AVISO SSH tracks only               Fcst 4
%        #6 - 3D T/S 1/30 NEMO
%        #7 - AVISO + UGOS PIES (T/S profiles)    Fcst 1
%        #8 - AVISO + extended PIES               Fcst 3
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
f_mhd_tser = 1; % plot MHD time series 


pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
btx = 'anls_fcst_mhd_LCLCEcntr.m';

iFcst  = 6;  % Forecast group to analyze, 1 - uses hindcast#7, 2 - #8, 3 - #2, 4 - #3
Nruns  = 1;  % # of runs in 1 forecast series
Ntimes = 1;  % time windows when fcst is initialized: 1 - May/June 2011, 2- Jan/Feb 2012

if iFcst<10
  pthmat = pthmat1;
else
  pthmat = pthmat2;
end

ixx=0;
ymx = 0;
irun1=1;
irun2=Nruns;
for ifcst=iFcst:iFcst
  switch (ifcst)
   case(1)
    nhnd=7;
   case(2)
    nhnd=8;
   case(3)
    nhnd=2;
    Ntimes=1;
    irun1=2;
    irun2=2;
   case(4);
    nhnd=3;
   case(6)
    nhnd=6;
   case(10)
    nhnd=10;   
  end
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



CLR = [0.4 0.8 0.2;...
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
clear ixx
for nhnd=6:6
  for itime=1:Ntimes  % 2011 and 2012 time windows
    for irun=1:1
      ifn=ifn+1;
      itot = (nhnd-7+itime-1)*Nruns+irun;

      if f_mhd_tser==1
        sub_plot_mhd_tser(MHD,CLR,irun,ifn,ipst,itot);
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
nll = length(MHD);
cMHD=[];
cMHDp=[];
icc=0;
for ill=1:nll
  if MHD(ill).iFcst ~= iFcst; continue; end;
  icc=icc+1;
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



figure(40); clf;

for itime=1:2 % summary for 2 time period of forecasts
  if itime==1
    axes('Position',[0.09 0.6 0.7 0.32]);
    j1=1;
    j2=7;
  else
    axes('Position',[0.09 0.1 0.7 0.32]);
    j1=8;
    j2=14;
  end
		hold on;
		for jm=1:3
				cfc = cMHD(j1:j2,jm);
				cps = cMHDp(j1:j2,jm);

				mFc = mean(cfc);
				Fc1 = min(cfc);
				Fc2 = max(cfc);

				mPs = mean(cps);
				Ps1 = min(cps);
				Ps2 = max(cps);

				clr1=[0.5 0.2 1];
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
				clr2=[0.7 0.7 0.7];
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
  ds1=MHD(j1).Name;
  ds2=MHD(j2).Name;
  dstr1=MHD(j1).Date_str;
  dstr2=MHD(j2).Date_str;
		stl=sprintf('MHD HYCOM/NEMO LCLCE %s/%s, Ini:%s-%s',ds1,ds2,dstr1,dstr2);
  title(stl);
end
axes('Position',[0.8 0.7 0.15 0.15]);
hold on
patch([1 1. 1.2 1.2],[1 1.2 1.2 1],clr1,'Edgecolor','none');
patch([1 1. 1.2 1.2],[0.7 0.9 0.9 0.7],clr2,'Edgecolor','none');
text(1.3,1.1,'Fcst','Fontsize',12);
text(1.3,0.8,'Persist','Fontsize',12);
set(gca,'xlim',[1 2],...
        'ylim',[0.6 1.2],...
        'visible','off');

bottom_text(btx,'pwd',1);		

	






