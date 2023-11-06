%  HYCOM forecasts
%
% Use T at depth to identidy LC - T and Z should be
%   same as for nemo: extr_lc_temp_nemo.m
%
% First, need to interpolate HYCOM T fields onto Z0
% interpolated in hycom_TSIS/distance_metric/interp_fcst_hycom2z.m
%
% Extract and save LC contours from HYCOM-TSIS and
%
% Compare LC contours from the HYCOM_TSIS forecasts
% using HYCOM hindcasts at day=1
% assimilating NEMO 1/100 free running fields
% and NEMO LC contours
%
% iFcst =   H/cast #2 - Full 2D SSH
%                 #3 - AVISO SSH tracks only
%                 #6 - 30th pnt T/S profiles GoM           not finished
%                 #7 - AVISO + UGOS PIES (T/S profiles)
%                 #8 - AVISO + extended PIES
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 1; % save mat; =2 - load saved and finish missing dates
% Set flags for extracting experiments:

iFcst = 8;
T0 = 2.5;   % T contour - anomaly demeaned T
Z0 = -200;  % depth


% Hindcasts:
pthd1  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
pthd2  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
pthd12 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_obs/';  % all satellites included
pthd3  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/2011_GLfreerun/'; % free run
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';


FCST = sub_fcst_info(iFcst);

%Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
HH(isnan(HH))=100;


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



% 1 forecast group at a time
Nruns=2;
nruns=Nruns;
Nhind = FCST.Nhind;  % initial cond from the  hindcast #
hnd_name = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
irun1 = FCST.run1;
irun2 = FCST.run2;
ntime=FCST.ntime; % 2 time windows for forecasts
for itime=1:ntime
  for irun=irun1:irun2
    RUN = FCST.TIME0(itime).RUN(irun);
    pthd1 = RUN.pthbin;
    TM    = RUN.TM;
    YDAY  = RUN.jday;
    nrc   = length(TM);
    DV    = datevec(TM);

    nmexp = sprintf('fcst%2.2i-%2.2i%2.2i',Nhind,itime,irun);
    ft2z = sprintf('%shycom_t2Z%4.4i_fcst%2.2i-%2.2i%2.2i.mat',...
                    pthmat,abs(Z0),Nhind,itime,irun);
    fprintf('\n\n %s %s\n',nmexp,ft2z);
    fprintf('Hindcast Name, iinitial fields: %s\n',hnd_name);
    fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',Nhind,itime,irun);
    fprintf(' Input data: %s\n',pthd1);
    fprintf(' Mat file: %s\n',ft2z);

    load(ft2z);
    TM_fcst = TZH.TM;
    nrc=length(TM_fcst);

% Output:
    fmatout = sprintf('%sHYCOM_dT%2.2i_Z%3.3i_contour_%s.mat',...
                    pthmat,round(T0*10),abs(Z0),nmexp);

    clear ii LCXY 
    cntr=0;
    for irc=1:nrc
      dnmb=TM_fcst(irc);
			dv   = datevec(dnmb);
			yr   = dv(1);
			mo   = dv(2);
			dm   = dv(3);

			tic;

			fprintf('Processing %s\n',datestr(dnmb));
			FF = squeeze(TZH.Tz(irc,:,:));

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

			cntr=cntr+1;
			LCXY.TM(cntr)    = dnmb;
			LCXY.XY(cntr).X  = XX;
			LCXY.XY(cntr).Y  = YY;

			fprintf('Processed 1 rec, %6.4f min\n\n',toc/60);
%keyboard
			if mod(irc,30)==0 & f_mat>0
					fprintf('Saving %s\n',fmatout);
					save(fmatout,'LCXY');
			end

		end;  % irun

		fprintf('Finished, Saving %s\n',fmatout);
		save(fmatout,'LCXY');

  end  % itime
end


