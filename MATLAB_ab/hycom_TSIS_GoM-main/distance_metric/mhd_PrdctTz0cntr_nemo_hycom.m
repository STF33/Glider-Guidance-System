% Predictability experiments
%
% Calculate MHD for the LC and LCE contours 
% in the Eastern GoM
% East of 90W
%
% for the HYCOM hindcasts and NEMO
% (1) interpolated in hycom_TSIS/distance_metric/interp_fcstPrdct_hycom2z.m
% (2) extracted in extr_lc_tempPrdct_hycom.m
%
% iFcst =   H/cast #2 - Full 2D SSH
%                 #3 - AVISO SSH tracks only
%                 #6 - T/S GoM at every 30th pnt on NEMO grid - only 2011
%                 #7 - AVISO + UGOS PIES (T/S profiles)
%                 #8 - AVISO + extended PIES
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

f_mat=1;


iFcst = 16; % #10, 11, ..., 16
T0 = 2.5;  % dt contour
Z0 = -200;
lnW = -90; % cutoff longitude

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

FCST  = sub_fcstPrdct_info(iFcst);

%
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



% YC point:
x0  = -85;
y0  =  22;

%if f_mhd==1
% Array with NEMO LC
fTNM = sprintf('%sNEMO_dT%2.2i_Z%3.3i_contour.mat',pthmat1,round(T0*10),abs(Z0));
fprintf('Loading %s\n',fTNM);
load(fTNM);
LCN=LCXY;
TMN=LCN.TM;


%TMN = LCN(1).TM;  % Nemo
nrc = length(LCN(1).XY);

% 1 forecast group at a time
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

    nmexp   = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
    fcntrTz = sprintf('%shycomPrdct_dT%2.2i_Z%3.3i_contour_fcst%2.2i-%2.2i%2.2i.mat',...
                    pthmat,round(T0*10),abs(Z0),iFcst,itime,irun);
    fprintf('\n\n %s %s\n',nmexp,fcntrTz);
    fprintf('Hindcast Name, iinitial fields: %s\n',hnd_name);
    fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',iFcst,itime,irun);
    fprintf(' Input data: %s\n',pthd1);
    fprintf(' Mat file: %s\n',fcntrTz);

    load(fcntrTz);


    TM_fcst = LCXY.TM;
    LCH=LCXY;   % HYCOM LC contours
    nrc=length(TM_fcst);
    clear Xh0 Yh0 MHD

    for irc=1:nrc
      dnmb=TM_fcst(irc);
      inemo=find(TMN==dnmb);

      if isempty(inemo),
        error('NEMO date not found %s\n',datestr(dnmb));
      end

						fprintf('Calculating MHD %s\n',nmexp);
						if mod(irc,100)==0,
								fprintf('  %6.2f%% done ...\n',irc/nrc*100);
						end

						dm1= TMN(inemo);
						ihc= irc;

						Xn = LCN.XY(inemo).X; % NEMO LC contour
						Yn = LCN.XY(inemo).Y;

% Clean contours near Cuba:
						xc1=-86.8;
						xc2=-84.07;
						xc3=xc2;
						xc4=-81.1;
						yc1=22.8;
						yc2=yc1;
						yc3=25.1;
						yc4=yc3;
						I=find(Yn<=yc1 & Xn>xc1);
						Xn(I)=nan;
						Yn(I)=nan;
						I=find(Yn<=yc3 & Xn>xc3);
						Xn(I)=nan;
						Yn(I)=nan;
						I=find(Xn<lnW);
						Xn(I)=nan;
						Yn(I)=nan;
						Inn = find(~isnan(Xn));
						Xn=Xn(Inn);
						Yn=Yn(Inn);

						Xh = LCH.XY(irc).X;
						Yh = LCH.XY(irc).Y;
%
      if isempty(Xh)
        fprintf('HYCOM contour is empty \n');
        keyboard;
      end

						I=find(Yh<=yc1 & Xh>xc1);
						Xh(I)=nan;
						Yh(I)=nan;
						I=find(Yh<=yc3 & Xh>xc3);
						Xh(I)=nan;
						Yh(I)=nan;
						I=find(Xh<lnW);
						Xh(I)=nan;
						Yh(I)=nan;
						Inn = find(~isnan(Xh));
						Xh=Xh(Inn);
						Yh=Yh(Inn);

% Keep persistence
      if (irc==1)
        Xh0=Xn;
        Yh0=Yn;
      end;

						P = [Xn,Yn];
						Q = [Xh,Yh];
						mhd1 = modified_hausdorff_distance(P,Q,'geo');

						P = [Xn,Yn];
						Q = [Xh0,Yh0];
						mhd0 = modified_hausdorff_distance(P,Q,'geo');

  				MHD(irc,1)  = mhd1;
  				MHD(irc,2)  = mhd0; % persistence
    end  % for irc
%keyboard				
				fprintf('  Fcst  min/max MHD=%8.2f  %8.2f\n',min(MHD(:,1)),max(MHD(:,1)));
    fprintf('  Perst min/max MHD=%8.2f  %8.2f\n',min(MHD(:,2)),max(MHD(:,2)));
				fmat1 = sprintf('%sMHD_dT%2.2i_Z%3.3i_hycom_Prdct%s.mat',...
																	pthmat,round(T0*10),abs(Z0),nmexp);
				fprintf('Saving %s\n',fmat1);
				save(fmat1,'MHD','TM_fcst');

		end    % for irun
	
end;
%end

fprintf('All done\n');


