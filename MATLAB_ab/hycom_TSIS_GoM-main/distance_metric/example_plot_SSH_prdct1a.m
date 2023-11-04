% HYCOM-TSIS forecast experiments
% epxeriments #: 10 - 16 corresponding to 
% start time in May/June 2011 (itime = 1) or Jan/Feb (itime = 2)
%
% There are 5 "runs" within each experiment
% irun = 1 - control run (the one that is used for comparison)
% Experiment # 10 itime = 1 irun = 1 starts on May 1, 2011 (NEMO initial conditions)
% run for 100 days
% Every other experiment is 7 day later 
%
% 4 perturbation experiments: use time t0 +/- dt=[1, 2] days
% where t0 is day #8 in the control run
%
% Saved ssh fields are for the GoM domain (not the whole HYCOM-TSIS domain)
% contact: FSU COAPS, D. Dukhovskoy
%
% Specify forecast, itime, irun:
close all
clear

% Forecast - initialized from hindcasts
% iFcst = 
% #10 - New set of forecasts initialized from 3D interpolated NEMO+GLORYS onto
%      HYCOM, see: anls_mtlb_utils/hycom_TSIS/interp_nemo
% #10 - May 1, 2011 and Jan 1, 2012
% #11 - May 8, 2011 and Jan 8, 2012
% etc
iFcst = 10;   % =10, ..., 16 - 7 runs 1 week apart
irun  = 1;    % irun1 - control f/cast initialized from NEMO - do not change irun
itime = 1;    % =1 May 2011, 2 = Jan 2012
fday  = 1;   % f/cast day
Z0    = -200; % cutoff depth for calculating mean SSH 

%pthmat = '/Net/kronos/ddmitry/hycom/TSIS/data2ncsu/';
pthmat = '/Net/ftp_pub/ddmitry/UGOS1_rmse/';

nmexp = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
fmatout = sprintf('%shycom_SSHGOM_Prdctb1a_fcst%2.2i-%2.2i%2.2i.mat',pthmat,iFcst,itime,irun);

load(fmatout);

GOM = [1   465
			 1     1
		 266     1
		 303    78
		 364   110
		 447   138
		 479   162
		 561   167
		 557   376
		 497   465];



TM   = SSH.Time;
dv   = datevec(TM(fday));
dnmb = TM(fday);
Iocn = SSH.Iocn;
LON  = SSH.LON;
LAT  = SSH.LAT;
HH   = SSH.HH;
dmm  = SSH.ssh(:,fday);


[mh,nh]=size(HH);
[XM,YM]=meshgrid([1:nh],[1:mh]);
INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
clear XM YM


ssh       = HH*nan;
ssh(Iocn) = dmm;

% Demean SSH for finding the LC contour:
ssh(INH==0)=nan;
I = find(HH<=Z0 & ~isnan(ssh));
ssh_mn = mean(ssh(I));
ssh = ssh-ssh_mn;

figure(1); clf;
pcolor(LON,LAT,ssh); shading flat; 
colorbar
caxis([-0.5 0.5]);
title(sprintf('%s, demeaned SSH (m), %s, meanSSH=%6.4fm',nmexp,datestr(dnmb),ssh_mn));



