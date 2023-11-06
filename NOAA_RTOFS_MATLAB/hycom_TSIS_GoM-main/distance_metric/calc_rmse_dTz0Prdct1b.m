% RMSE dT at depth Z0  
% experiments 1b - compare hycom perturbed f/casts 
% to the hycom control f/cast
% contolr run is done for 1a experiments (initialized from 
%  interpolated NEMO+GLORYS fields run for 100 days)
%
% f/casts initialized from interpolated NEMO+GLORYS
%
% Calculate RMSE between NEMO and 
% HYCOM analysis SSH fields
% Predictability experiments: iFcst 10, 11, ..., 16:
% Forecasts: use interpolated fields as initial fields
%            within each f/cast - perorm additional f/cast runs
%            initialized from day t0=7 days of HYCOM main f/cast
%            run1 = main f/cast, 
%            run2 = -2day shift used as IC at day t0
%            run3 = -1 day shift, etc.
%            run4 = +1 days
%            run5 = +2 days
%


addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

%           
iFcst = 16; % iFcst = hindcast # used for initial fields
irun1=2;   % irun1 - control f/cast, 2,.., 5 - shifted +/- 2, 1 perturbation f/casts
irun2=5;
itime1=1;   % =1 May 2011, 2 = Jan 2012
itime2=2;   % 

f_mat = 1; % =1 save and overide existing mat; =2 - load saved and finish missing dates

if f_mat==0
  fprintf('Output is not saved !!! \n\n');
  fprintf('Check f_mat=%i\n',f_mat);
end

Z0 = -200; % discard close to coastline points

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';
btx      = 'calc_rmse_dTz0Prdct1b.m';

% Info about All hindcast experiments
load('hycom_tsis_expts.mat');

%ts times:
FCST = sub_fcstPrdct_info(iFcst);

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


% Indices for subsampled GoM domain 
[XX,YY,HHs,INH,Iocn,xt1,xt2,yt1,yt2,its1,its2,jts1,jts2] = ...
   sub_subsample_HYCOM2GoM(HH,LON,LAT,Z0);
Nocn=length(Iocn);

% 1 forecast group at a time
nmexp = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
for itime=itime1:itime2
  for irun=irun1:irun2    % forecast runs
    clear TZH
    ft2z = sprintf('%shycom_t2Z%4.4i_fcstPrdct1b%2.2i-%2.2i%2.2i.mat',...
                    pthmat,abs(Z0),iFcst,itime,irun);
    fprintf('Loading T2Z fields %s\n',ft2z);
    load(ft2z);

    TM   = TZH.TM;
    nrc  = length(TM);
    DV   = datevec(TM);
    cntr = 0;
    LONH = TZH.LON;
    LATH = TZH.LAT;
				HHx  = TZH.HH; % same as HHs

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
    frmseout = sprintf('%sRMSE_Prdct1b_dT%4.4i_%s.mat',pthout,abs(Z0),fcstname);

    fprintf('\n\n %s %s\n',nmexp,frmseout);
				fprintf(' Input data: %s\n',ft2z);


% Get persistence - day 1 from NEMO:
    dnmb  = TM(1);
    irun0 = 1;   % control run, irun = 1
    dmm   = sub_getTz0_hycom_control(dnmb,iFcst,itime,irun0,abs(Z0));
    [thmP,tMhP]   = sub_demeanTz_hycom(dmm);
    IocnTz0 = find(~isnan(thmP)); % ocean points inside GoM and 200m
    [mmh,nnh] = size(HHx);

    clear RMSERR
    RMSERR.Name       = EXPT(iFcst).Name;
    RMSERR.Name_short = EXPT(iFcst).Name_short;
    RMSERR.Indx_gom   = IocnTz0;
    RMSERR.HYCOM_gom  = [mmh,nnh];
    RMSERR.HH_gom     = HHx;    

%keyboard
% Continue from the last saved record
% skip finished mat files
    last_rec=0;
    if f_mat==2
      fprintf('Continue from the last unfinished mat file\n');
      if exist(frmseout,'file');
        fprintf('Loading %s\n',frmseout);
        load(frmseout);
        [a1,a2]=size(RMSERR.ERR_squared);
        last_rec=a2;
%        keyboard
      end
    end
        

%				cntr=0;
				for ii=1:nrc
      if f_mat==2
        if ii<=last_rec 
          fprintf('%s Record exist skipping %i/%i/%i\n',fcstname,DV(ii,1:3));
          continue
        else
          fprintf('Continue from ii=%i %s %i/%i/%i\n\n',ii,fcstname,DV(ii,1:3));
          f_mat=1;
          last_rec=0;
        end
      end
%      keyboard
						tic;

						yr  = DV(ii,1);
						mo  = DV(ii,2);
						dm  = DV(ii,3);
						dnmb= TM(ii);

% HYCOM control run:
      FF0 = sub_getTz0_hycom_control(dnmb,iFcst,itime,irun0,abs(Z0));
      [thm0,tMh0] = sub_demeanTz_hycom(FF0);


% HYCOM field perturbed field:
      FF = squeeze(TZH.Tz(ii,:,:));
      [thm,tMh] = sub_demeanTz_hycom(FF);

      f_pltT = 0;
      if f_pltT==1
        figure(20); clf;
        hold on;
        pcolor(LONH,LATH,thm0); shading flat;
        contour(LONH,LATH,thm0,[2.5 2.5],'k-');
        axis('equal');

        colorbar
        caxis([-4 4]);
        stl = sprintf('HYCOM Control, dT z=%i <T>=%5.2f %i/%i/%i',abs(Z0),tMh0,dm,mo,yr);
        title(stl);
        bottom_text(btx,'pwd',1);

        figure(21); clf;
        hold on;
        pcolor(LONH,LATH,thm); shading flat;
        contour(LONH,LATH,thm,[2.5 2.5],'k-');
        axis('equal');

        colorbar
        caxis([-4 4]);
        stl = sprintf('HYCOM %s, dT z=%i <T>=%5.2f %i/%i/%i',fcstname,abs(Z0),tMh,dm,mo,yr);
        title(stl);
        bottom_text(btx,'pwd',1);
        keyboard
      end       
 

      RMSE=[];
      RMSE  = (thm(IocnTz0)-thm0(IocnTz0)).^2;  % f/cast - control run
      RMSEP = (thmP(IocnTz0)-thm0(IocnTz0)).^2;  % f/cast - persistence
      cnn = length(IocnTz0);
      rms=sqrt(nansum(RMSE)/cnn);
      rmsp=sqrt(nansum(RMSEP)/cnn);

% Spatial map of RMSE;
      f_chck=0;
      if f_chck==1
        ERR=zeros(mmh,nnh)*nan;
        ERR(IocnTz0)=sqrt(RMSE);
        figure(11); clf;
        pcolor(ERR); shading flat;
        set(gca,'xlim',[1 nnh],...
                'ylim',[1 mmh]);
        colorbar

        keyboard
      end

      RMSE=RMSE(:);
      RMSERR.ERR_squared(:,ii)=RMSE;
      RMSERR.ERRprst_squared(:,ii)=RMSEP;  % persistence squared err
      RMSERR.MeanGoM_T_hycom(ii,1) = tMh; 
      RMSERR.MeanGoM_T_perst(ii,1) = tMhP; % persistence mean T

      fprintf('1 record %6.3f min, RMSE=%6.4g, PerstRMSE=%6.4g, MnT=%4.1f MnT_perst=%4.1f \n',...
              toc/60,rms,rmsp,tMh,tMhP);

      if mod(ii,15)==0 & f_mat==1  % f_mat=2 should change to f_mat=1
        fprintf('Saving %s\n',frmseout);
        save(frmseout,'RMSERR');
      end
    end

    if f_mat==1
      fprintf('Saving %s\n',frmseout);
      save(frmseout,'RMSERR');
    end

  end % irun loop - 
end

