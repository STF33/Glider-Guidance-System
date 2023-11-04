function FCST = sub_fcst_info(iFcst);
% Provide information about forecast iFcst
% iFcst = the hindcast number used for intialization
% All experiments <10 - OSSEs experiments with initial fields from IAS HYCOM-TSIS assimilated
%     hindcast runs
% Experiments 10 - new predictability experiments with initial fields from 
%             interpolated NEMO+GLORYS fields onto HYCOM
%
pthfrcst  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';
%pthfrcst2 = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_predictability/';
pthfrcst2 = '/nexsan/people/ddmitry/hycom/TSIS/IASx0.03/output_predictability/';

FCST = struct;
cc=0;
dd=7;  % days shift in forecast start date
Nruns=7; % # of forecast runs in each time window
ntime=2; % 2 time windows/periods (LC state) for forecasts 2011/2012
irun1=1;
irun2=Nruns;
if iFcst>=10; 
  error('For Predictability runs iFcst>=10 use sub_fcstPrdct_info.m');
end; % not all f/casts finished

imm=0;

% All hindcast experiments
load('hycom_tsis_expts.mat');

% Forecast - initialized from hindcasts
% iFcst = 
% #2 - full 2D SSH
% #3 - AVISO ssh swaths
% #6 - T/S GoM at every 30th pnt on NEMO grid
% #7 - AVISO Swath UGOS PIES
% #8 - AVISO Swath extended UGOS
%
% Predictability experiments - start from interpolated NEMO+GLORYS
%  - usei sub_fcstPrdct_info.m

nhind=iFcst;
ntime=2;
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
% Different forecast time windows May.June 2011 and Jan/Feb 2012
for idd=1:ntime
		switch(idd)
			case(1)
				day0 = datenum(2011,05,01);
			case(2)
				day0 = datenum(2012,01,01);
		end

		dnmb=day0-dd;

		for irun=1:Nruns  % forecast runs, 90 days each
				dnmb=dnmb+dd;
				dv=datevec(dnmb);

				drnm = sprintf('fcst%2.2i_%2.2i%2.2i%4.4i',nhind,dv(3),dv(2),dv(1));
				pthbin = sprintf('%s%s/',pthfrcst,drnm);
				fprintf('Forecast path: %s\n',pthbin);

				FCST.TIME0(idd).RUN(irun).pthbin = pthbin;
				id1 = dnmb-datenum(dv(1),1,1)+1;
				icc=0;
				FCST.TIME0(idd).RUN(irun).jday = [id1:id1+90];
				FCST.TIME0(idd).RUN(irun).TM   = [dnmb:dnmb+90];
		end
end



return
