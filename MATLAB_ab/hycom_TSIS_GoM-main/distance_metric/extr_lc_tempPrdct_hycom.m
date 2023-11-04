% Predictability experiments
% interpolated NEMO+GLORYS onto HYCOM grid
% used as initial fields 
%
% Analysis of LC, steps:
% First, need to interpolate HYCOM T fields onto Z0
%        interp_hycom2z.m
%
% 1) extract LC contour extr_lc_tempPrdct_hycom.m
% 2) same for nemo LC contour:  extr_lc_temp_nemo.m if needed
% 3) calculate MHD: distance_metric/mhd_osse_hindcasts_hycom_nemoV1.m
% 4) Plot results: distance_metrics/
%
% Use T at depth to identidy LC - T and Z should be
%   same as for nemo: extr_lc_temp_nemo.m
%
% Compare LC contours from the HYCOM_TSIS hindcast
% assimilating NEMO 1/100 free running fields
% and NEMO LC contours
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 1; % save mat; =2 - load saved and finish missing dates
% Forecast - initialized from hindcasts
% iFcst = 
% #10 - May 1, 2011 (Time=1) and Jan 1, 2012 (Time=2)
% #11 - May 8, 2011 and Jan 8, 2012
% etc
iFcst=16; 

T0 = 2.5;   % T contour - anomaly demeaned T
Z0 = -200;  % depth


% Hindcasts:
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';


FCST  = sub_fcstPrdct_info(iFcst);
ntime = FCST.ntime; % 2 time windows for forecasts
irun1 = FCST.run1;
irun2 = FCST.run2;

%
% HYCOM:
rg=9806;  % convert pressure to depth, m
huge=1e20;
btx = 'extr_lc_tempPrdct_hycom.m';


% GoM region HYCOM:
% Note - indices for subsampled GoMregion
% from HYCOM TSIS IAS
GOM=[362   108
   288    50
   159     4
     9     5
     3   141
     8   401
   330   407
   398   405
   511   305
   520   210
   507   148
   429   118];



% Read in HYCOM ssh from the forecast runs
nfcst = length(FCST);
fcst_name = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
irun1 = FCST.run1;
irun2 = FCST.run2;
ntime = FCST.ntime;
for itime=1:ntime
  for irun=irun1:irun2
    clear TZH
    ft2z = sprintf('%shycom_t2Z%4.4i_fcstPrdct%2.2i-%2.2i%2.2i.mat',...
                    pthmat,abs(Z0),iFcst,itime,irun);
    fprintf('Loading T2Z fields %s\n',ft2z);
    load(ft2z);

% Output:
    fmatout = sprintf('%shycomPrdct_dT%2.2i_Z%3.3i_contour_fcst%2.2i-%2.2i%2.2i.mat',...
                    pthmat,round(T0*10),abs(Z0),iFcst,itime,irun);


				TM  = TZH.TM;
				nrc = length(TM);
				cntr=0;

				for ii=1:nrc
						dnmb = TM(ii);

%						if dnmb<1000; continue; end % indexing bug in interp_hycom2z.m started from unfinished mat
						dv   = datevec(dnmb);
						yr   = dv(1);
						mo   = dv(2);
						dm   = dv(3);

						tic;

						fprintf('Processing %s\n',datestr(dnmb));
						FF = squeeze(TZH.Tz(ii,:,:));

						if ~exist('INH','var');
								HHx = TZH.HH;
								[mm,nn]=size(HHx);
								[XM,YM]=meshgrid([1:nn],[1:mm]);
								INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
								clear XM YM;

								LONH = TZH.LON;
								LATH = TZH.LAT;
						end

				% Subtract spatial mean T
						dmm=FF;
						dmm(INH==0)=nan;
						tM=nanmean(nanmean(dmm));
						tnm = FF-tM;
						tnm(INH==0)=nan;

						dmm=tnm;
						dmm(INH==0)=nan;

						Bisol=T0;
						npnt0=200; % min # of points for the contour
						[XX,YY] = get_contour(LONH,LATH,dmm,Bisol,npnt0);
%keyboard
						cntr=cntr+1;
						LCXY.TM(cntr)    = dnmb;
						LCXY.XY(cntr).X  = XX;
						LCXY.XY(cntr).Y  = YY;

						fprintf('Processed 1 rec, %6.4f min\n\n',toc/60);

						if mod(ii,30)==0 & f_mat>0
								fprintf('Saving %s\n',fmatout);
								save(fmatout,'LCXY');
						end

				end;

    fprintf('Finished, Saving %s\n',fmatout);
    save(fmatout,'LCXY');

  end
end


