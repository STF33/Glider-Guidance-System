% Subsample SSH fields to a smaller domain
% Experiment Predictability 1b
% For NCSU group

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 1;

%Z0 = -10; % discard close to coastline points
%Z0 = -200;

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/data2ncsu/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';


% Forecast - initialized from hindcasts
% iFcst = 
% #10 - New set of forecasts initialized from 3D interpolated NEMO+GLORYS onto
%      HYCOM, see: anls_mtlb_utils/hycom_TSIS/interp_nemo
% #10 - May 1, 2011 and Jan 1, 2012
% #11 - May 8, 2011 and Jan 8, 2012
% etc
%iFcst=10;  % =10, ..., 16 - 7 runs 1 week apart
ifcst1=10;
ifcst2=16;
irun1=1;   % irun1 - control f/cast, 2,.., 5 - shifted +/- 2, 1 perturbation f/casts
irun2=1;
itime1=1;   % =1 May 2011, 2 = Jan 2012
itime2=2;   % 

f_mat = 1; % >0 save mat;
% Hindcasts:
%  pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_predictability/';
pthfcst = '/nexsan/people/ddmitry/hycom/TSIS/IASx0.03/output_predictability';
%pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/data2ncsu/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';

% All hindcast experiments
load('hycom_tsis_expts.mat');

% iFcst=10,..., 16
%ntime=FCST.ntime; % 2 time windows for forecasts
%irun1=FCST.run1;
%irun2=FCST.run2;

% HYCOM:
rg=9806;  % convert pressure to depth, m
huge=1e20;

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
Iocn = find(INH==1 & HH<0);
clear XM YM


% Subsample to a smaller domain:
xnas1 = min(LON(Iocn));
xnas2 = max(LON(Iocn));
ynas1 = min(LAT(Iocn));
ynas2 = max(LAT(Iocn));
[J,I]=find(LON<xnas1);
its1=max(I);
[J,I]=find(LON>xnas2);
its2=min(I);
[J,I]=find(LAT<ynas1);
jts1=max(J);
[J,I]=find(LAT>ynas2);
jts2=min(J);

xt1=LON(jts1,its1);
yt1=LAT(jts1,its1);
xt2=LON(jts2,its2);
yt2=LAT(jts2,its2);

LON_gom=LON(jts1:jts2,its1:its2);
LAT_gom=LAT(jts1:jts2,its1:its2);
HH_gom=HH(jts1:jts2,its1:its2);
[mgom,ngom] = size(HH_gom);
IocnG = find(HH_gom<0);



for itime=itime1:itime2
  SSH = struct;
  for iFcst=ifcst1:ifcst2
    FCST = sub_fcstPrdct_info(iFcst);
    irun =  1;  % only 1 run in Prdct1a for eeach expt
		RUN0  = FCST.TIME0(itime).RUN(1); % control run
		RUN   = FCST.TIME0(itime).RUN(irun);
		pthd1 = RUN.pthbin;
		TM    = RUN.TM;
		YDAY  = RUN.jday;
		nrc   = length(TM);
		DV    = datevec(TM);

    nmexp = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
		fmatout = sprintf('%shycom_SSHGOM_Prdctb1a_fcst%2.2i-%2.2i%2.2i.mat',pthmat,iFcst,itime,irun);
%
		day0 = FCST.TIME0(itime).RUN(irun).prdct_time0;
		irc1=1;

    fprintf('\n\n %s %s\n',nmexp,fmatout);
    fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',iFcst,itime,irun);
    fprintf(' Input data: %s\n',pthd1);
    cntr=0;
    irc=0;
    for ii=1:nrc
      if ii<irc1, continue; end;
			yr  = DV(ii,1);
			mo  = DV(ii,2);
			dm  = DV(ii,3);
			dnmb= TM(ii);
			iday= YDAY(ii);

% Day 1 - initial field - take from the original f/cast if needed for perturbed cases
			fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
			finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
			fin=fina;

			ie = exist(fin,'file');

			if ~ie
				fprintf('  ERR ** Missing forecast: %s\n',fin);
				keyboard
				continue;
			end

			fprintf('Reading %s\n',fina);
			fld = 'srfhgt';
			[F,nn,mm,ll] = read_hycom(fina,finb,fld);
			F(F>huge)=nan;
			ssh=squeeze(F)./(1e-3*rg);  % ssh m
% Subsample
      dmm = ssh(jts1:jts2,its1:its2);
      ssh_gom = dmm(IocnG);

      irc=irc+1;
      SSH.ssh(:,irc)=ssh_gom;
    end  % irun
% save output
		SSH.Time = TM;
		SSH.Iocn = IocnG;
		SSH.HH   = HH_gom;
		SSH.LON  = LON_gom;
		SSH.LAT  = LAT_gom;
%
%keyboard
    if f_mat==1
      fprintf('Saving %s\n\n',fmatout);
      save(fmatout,'SSH');
    end

  end
end


