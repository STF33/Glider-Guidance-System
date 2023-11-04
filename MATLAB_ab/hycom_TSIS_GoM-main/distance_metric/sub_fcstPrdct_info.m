% Provide information about forecast iFcst
% for Predictability experiments:
% Predictability experiments - start from interpolated NEMO+GLORYS
% Initialize 100-day f/casts same dates / time periods as in the hindcast/forecasts experiments
% i.e. TPeriod 1: May 1,  2011 -  #10 f/cast
%                 May 8,  2011 -  #11 
%                 May 15, 2011 -  #12
%                 ...
% 
%      TPeriod 2: Jan 2, 2012  -  #10
%                 Jan 9, 2012  -  #11
%   etc
%
% run1 - original run (100 days), use t0 = day 8 from run1 as unperturbed run
% 
%  Each forecast has 4 perturbation runs initialized from the forecast field 
%       at day t0 = 7 +/- 1,2 days
%
% e.g.: fcst10, TPeriod1: run 2 = use fields from May 06 as May 08
%                         and integrate for 90 days
%                     run3: May 07 --> May 08 (consider May 07 as May 08)
%                     already exists: May 08 --> May 08 - control run
%                     run4: May 09 --> May 08
%                     run 5: May10 --> May 08
%
function FCST = sub_fcstPrdct_info(iFcst);

pthfrcst  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';
%pthfrcst2 = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_predictability/';
pthfrcst2 = '/nexsan/people/ddmitry/hycom/TSIS/IASx0.03/output_predictability/';

FCST = struct;
cc=0;
dd=7;  % days shift in main forecast start date
Nruns=5; % # of forecast perturbation runs in each main f/cast and Time Period
ntime=2; % 2 time windows/periods (LC state) for forecasts 2011/2012
irun1=1;
irun2=Nruns;


imm=0;

% All hindcast experiments
load('hycom_tsis_expts.mat');

% Predictability experiments - start from interpolated NEMO+GLORYS
%
% #10 - May 1, 2011 / Jan 1, 2012 2 Time periods + 4 extra runs +/- 1,2 days
% #11 - May 8

nhind=iFcst;
%switch(iFcst)
%	case(2)
%		ntime=1; % only 2011 at the moment completed
% otherwise
%  ntime=2;
%end
%end
%imm=imm+1;
FCST.Nhind = nhind;
FCST.Hindcast_Name = EXPT(nhind).Name;
FCST.Hindcast_path = EXPT(nhind).path;
FCST.run1=irun1;
FCST.run2=irun2;
FCST.ntime=ntime;
%
% Different forecasts have different starting date
% for time windows May.June 2011 and Jan/Feb 2012
add_day = (iFcst-10)*dd;  % start day offset for 7 forecast groups (iFcst)
time0_prtrb = 7; % time 0 for the perturb experiments in the control run
dDay = [0; -2; -1; 1; 2]; % start day shift in the perturb experiments
for idd=1:ntime
		switch(idd)
			case(1)
				day0 = datenum(2011,05,01)+add_day; % day 1 of the control run
			case(2)
				day0 = datenum(2012,01,01)+add_day; % day 1 of the contolr run
		end

		for irun=1:Nruns  % forecast runs, 100 for control or 90 days for perturb
    dd = dDay(irun);

    if irun==1
	  	  dnmb=day0;
    else
	  	  dnmb=day0+time0_prtrb;
    end
				dv=datevec(dnmb);
    dnmbT0_true=day0+time0_prtrb+dd; % actual date0 of the f/cast at t=0
    dv_true = datevec(dnmbT0_true);
    

				drnm = sprintf('fcst%2.2i_%2.2i%2.2i%4.4i',nhind,dv(3),dv(2),dv(1));
    if irun>1
      drnm = sprintf('fcst%2.2ir%i_%2.2i%2.2i%4.4i',nhind,irun,dv(3),dv(2),dv(1));
    end
				pthbin = sprintf('%s%s/',pthfrcst2,drnm);
				fprintf('Forecast path: %s\n',pthbin);

				FCST.TIME0(idd).RUN(irun).pthbin = pthbin;
				id1 = dnmb-datenum(dv(1),1,1)+1;
				icc=0;
    if irun==1
      ndays=99;
    else
      ndays=89;
    end
				FCST.TIME0(idd).RUN(irun).jday = [id1:id1+ndays];
				FCST.TIME0(idd).RUN(irun).TM   = [dnmb:dnmb+ndays];
    FCST.TIME0(idd).RUN(irun).dayT0_true = dnmbT0_true;
    FCST.TIME0(idd).RUN(irun).prdct_time0 = day0+time0_prtrb;
		end
end



return
