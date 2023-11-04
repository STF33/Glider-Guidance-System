% Predictability experiments: 
% initial fields: NEMO+GLORYS interpolated onto HYCOM-TSIS
%
% Analysis of LC, steps:
% 1) extract LC contour from the forecast runs
% 2) same for nemo LC contour:  extr_lc_ssh_nemo.m if needed
% 3) calculate MHD: distance_metric/mhd_osse_hindcasts_hycom_nemoV1.m
% 4) Plot results: distance_metrics/
%
%  Forecasts: decision tree:
%  #10 - # 16, Time 1, 2 - main f/casts
%  then for each f.cast - +/-, 1, 2 days from day t0=7 of the main f/cast
%  'predturbation expts"
%      |
%      V
%  TimePeriod during which the f/casts are initialized (May-June 2011, Jan-Feb 2012)
%     |
%     V
%  Forecasts (90 day f/csts shifted 7 days - 7 fcsts)
%
%
% Extract and save LC contours from HYCOM control runs / forecasts 
% specify individually which run need to extract
% LC and LCE: saved 1 LCEs (the largest) if there is one
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

% Forecast - initialized from hindcasts
% iFcst = 2, 3, 6, 7, 8 - OSSE f/casts
% #10 - New set of forecasts initialized from 3D interpolated NEMO+GLORYS onto
%      HYCOM, see: anls_mtlb_utils/hycom_TSIS/interp_nemo
% #10 - May 1, 2011 and Jan 1, 2012
% #11 - May 8, 2011 and Jan 8, 2012
% etc
% 1 forecast at a time
iFcst=8; % #2, #3, #6, #7, #8 = match the hindcast numbers used for iniialization


f_mat = 0; % >0 save mat;

Bisol = 0.17;  % ssh contour
huge  = 2e20;  
rg    = 9806;


% Hindcasts:
%  pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_predictability/';
pthfcst = '/nexsan/people/ddmitry/hycom/TSIS/IASx0.03/output_predictability';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';

%fmatout = sprintf('%sLC_coord_osse_hycomV1.mat',pthmat);
%
% All hindcast experiments
load('hycom_tsis_expts.mat');



% Forecasts times:
FCST = struct;
cc=0;
dd=7;  % days shift in forecast start date
imm=0;

if iFcst<10
  FCST = sub_fcst_info(iFcst);
else
  FCST = sub_fcstPrdct_info(iFcst);
end
Nhnd = FCST.Nhind;  % initial cond from the  hindcast #
hnd_name = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
irun1 = FCST.run1;
irun2 = FCST.run2;
ntime=FCST.ntime; % 2 time windows for forecasts
NTime=ntime;

% All hindcast experiments
load('hycom_tsis_expts.mat');


% Array of forecast runs
icc=0;
for it=1:NTime
  for irun=irun1:irun2
    icc=icc+1;
    FRUNS(icc).Nhndcst = Nhnd;
    FRUNS(icc).TimePeriod = it;
    FRUNS(icc).Nfcst = irun;
    FRUNS(icc).matname = sprintf('hycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',Nhnd,it,irun);
    FRUNS(icc).Hind_name = EXPT(Nhnd).Name;
  end
end


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
clear XM YM


Mruns = length(FRUNS);
for ixx = 1:Mruns
  nmexp = FRUNS(ixx).Hind_name;
  Nhnd  = FRUNS(ixx).Nhndcst;
  flnm = FRUNS(ixx).matname;
  nmrun = flnm;
  fmatout = sprintf('%s%s',pthmat,flnm);
  pthHcst = EXPT(Nhnd).path;    % hindcast path with initial fields
  itime   = FRUNS(ixx).TimePeriod;
  irun    = FRUNS(ixx).Nfcst;

  RUN0  = FCST.TIME0(itime).RUN(1); % control run
  RUN   = FCST.TIME0(itime).RUN(irun);
  pthd1 = RUN.pthbin;
  TM    = RUN.TM;
  YDAY  = RUN.jday;
  nrc   = length(TM);
  DV    = datevec(TM);
      
  fmatout = sprintf('%shycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',pthmat,Nhnd,itime,irun);

  if exist(fmatout,'file') & f_mat>0
    fprintf('!!!  %s exist, skipping ...\n',fmatout);
    keyboard
%        continue
  end

  fprintf('\n\n %s %s\n',nmexp,fmatout);
  fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',Nhnd,itime,irun);
  fprintf(' Input data: %s\n',pthd1);

  clear LCXY
  LCXY.Name = nmexp;
  LCXY.Pthdata = pthd1;

  cntr=0;
  for ii=1:nrc
    tic;

    yr  = DV(ii,1);
    mo  = DV(ii,2);
    dm  = DV(ii,3);
    dnmb= TM(ii);
    iday= YDAY(ii);
  
    fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
    finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
% Day 1 (initial field) from the original hindcast
    if ii==1
      fprintf('Initial fields from hindcast\n');
      fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthHcst,yr,iday);
      finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthHcst,yr,iday);
    end 
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
  %
  % Subtract anomaly:
    dmm=ssh;
    dmm(INH==0)=nan;
    dmm(HH>-200)=nan;
    sshM=nanmean(nanmean(dmm));
    ssh=ssh-sshM;

  %
  % Derive LC contour:
  % 
    dmm=ssh;
    dmm(INH==0)=nan; 
    bsgn = -1;
    f_stop = 0;
    LCH1 = identify_LC(LON,LAT,dmm,Bisol,'err_stop',f_stop);
%    keyboard

    cntr=cntr+1;

%        if cntr==19;
%          keyboard;
%        end

  % HYCOM-TSIS
    LCXY.TM(cntr)    = dnmb;
    LCXY.XY(cntr).X  = LCH1(1).xx;
    LCXY.XY(cntr).Y  = LCH1(1).yy;
% Save LCEs as well:
    lcc=length(LCH1);  % LC + # of LCEs, first is LC
    LCE(1).NumbLCE(cntr) = lcc;
    LCE(1).TM(cntr)  = dnmb;
    LCE(1).XY(cntr).X=[];
    LCE(1).XY(cntr).Y=[];
    nlce = length(LCE);
    for ilce=2:nlce
      LCE(ilce).XY(cntr).X=[];  % LCE contours empty if no LCE's
      LCE(ilce).XY(cntr).Y=[];
    end

    if lcc>1
      for ilc=2:lcc
        LCE(ilc-1).XY(cntr).X = LCH1(ilc).xx;
        LCE(ilc-1).XY(cntr).Y = LCH1(ilc).yy;
      end
    end
    fprintf('cntr=%i, Processed 1 rec, %6.4f min\n\n',cntr,toc/60);

    btx = 'extr_lc_hycom_OSSEfcst.m';
    f_chck = 0;
    if ii>=20 & f_chck==1
      figure(10); clf;
      pcolor(LON,LAT,dmm); shading flat;
      hold on; 
      contour(LON,LAT,dmm,[Bisol Bisol],'k');

      axis('equal');
      set(gca,'xlim',[-98 -80],...
              'ylim',[17 31]);

      xlc = LCXY.XY(cntr).X;
      ylc = LCXY.XY(cntr).Y;
      dstr = datestr(dnmb);
      plot(xlc,ylc,'r.');

      nlce = length(LCE);
      for ilce=1:nlce
       xlce = LCE(ilce).XY(cntr).X;
       ylce = LCE(ilce).XY(cntr).Y;
       if isempty(xlce), continue; end;

       plot(xlce,ylce,'m.'); 
      end
      stl = sprintf('%s %s',nmrun,dstr); 
      title(stl);
      bottom_text(btx,'pwd',1);

      keyboard
    end

    if mod(ii,20)==0 & f_mat>0
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'LCXY','LCE');
    end
  end; % days in 1 run

  fprintf('Finished, Saving %s\n',fmatout);
  save(fmatout,'LCXY','LCE');
  % 
end



