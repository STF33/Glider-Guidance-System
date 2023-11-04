% RMSE dT at depth Z0  
% Predictability experiments
% HYCOM T should be interpolated onto depth layer z0 first
%
% Layer 1 in HYCOM 3 m thick
% NEMO: Layers 1 = -1.05 m, Layer 3 (bottom interface) = -3.26198
%
% Predictability experiments:
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
%
%      |
%      V
%  TimePeriod during which the f/casts are initialized (May-June 2011, Jan-Feb 2012)
%     |
%     V
%  Forecasts (90 day f/csts shifted 7 days - 7 fcsts)
%
%


addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

%  Predictability forecasts (initialized from NEMO+GLORYS interpolation):
%                 #10 - NEMO+GLORYS at day 1
%           
iFcst = 16; % iFcst = hindcast # used for initial fields
f_mat = 1; % =1 save and overide existing mat; =2 - load saved and finish missing dates

if f_mat==0
  fprintf('Output is not saved !!! \n\n');
  fprintf('Check f_mat=%i\n',f_mat);
end

Z0 = -200; % discard close to coastline points

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';
pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';

btx = 'calc_rmse_dTz0Prdct.m';

% NEMO 2 HYCOM indices:
%fmatout = sprintf('%snemo2hycom_indx.mat',pthmat);
%load(fmatout);


% Info about All hindcast experiments
load('hycom_tsis_expts.mat');

% Load NEMO 2 hycom indices:
fmatindx = sprintf('%snemo2hycom_indx.mat',pthmat);
fprintf('Loading %s\n',fmatindx);
load(fmatindx);

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


if m~=INDX.HYCOM_size(1) | n~=INDX.HYCOM_size(2)
  fprintf('Saved domain: %i %i, this domain: %i %i\n',...
   INDX.HYCOM_size(1:2), m, n)
  error('Saved NEMO indices are for different HYCOM domain');
end


IJhycom = INDX.HYCOM_ocean;  % saved HYCOM pnts with NEMO indices

[XX,YY,HHs,INH,Iocn,xt1,xt2,yt1,yt2] = sub_subsample_HYCOM2GoM(HH,LON,LAT,Z0);


% 
[LONN,LATN,INN] = sub_getNEMO_info;
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
fprintf('Loading NEMO grid %s\n',fgrd);
load(fgrd);
dZN = abs(diff(ZZN)); 
% Find depth for NEMO:
dZ = abs(abs(ZZN)-abs(Z0));
iz0 = find(dZ==min(dZ));


% 1 forecast group at a time
nmexp = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
irun1 = FCST.run1;
irun2 = FCST.run2;
ntime = FCST.ntime;
for itime=1:ntime
  for irun=irun1:irun2    % forecast runs
    clear TZH
    ft2z = sprintf('%shycom_t2Z%4.4i_fcstPrdct%2.2i-%2.2i%2.2i.mat',...
                    pthout,abs(Z0),iFcst,itime,irun);
    fprintf('Loading T2Z fields %s\n',ft2z);
    load(ft2z);

    TM  = TZH.TM;
    nrc = length(TM);
    DV  = datevec(TM);
    cntr=0;

    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',iFcst,itime,irun);
    frmseout = sprintf('%sRMSE_Prdct_dT%4.4i_%s.mat',pthout,abs(Z0),fcstname);

    fprintf('\n\n %s %s\n',nmexp,frmseout);
		fprintf(' Input data: %s\n',ft2z);

    clear RMSERR
    RMSERR.Name       = EXPT(iFcst).Name;
    RMSERR.Name_short = EXPT(iFcst).Name_short;
    RMSERR.Indx_hycom = Iocn;
    RMSERR.HYCOM_dim  = [mh,nh];

% Get persistence - day 1 from NEMO:
    dnmb=TM(1);
    tnm = sub_getTz0_nemo(dnmb,iz0);
    tMn=nanmean(nanmean(tnm));
    tnmP = tnm-tMn;
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
%
% Nemo T fields:
			tnm = sub_getTz0_nemo(dnmb,iz0);
			tMn=nanmean(nanmean(tnm));
			tnm = tnm-tMn;
% HYCOM field:
			FF = squeeze(TZH.Tz(ii,:,:));

			if ~exist('joff','var') | isempty(joff);
				HHx = TZH.HH;
				[mm,nn]=size(HHx);
				[XM,YM]=meshgrid([1:nn],[1:mm]);
				LONH = TZH.LON;
				LATH = TZH.LAT;
	%
				Nocn=length(Iocn);
%
% Find offset for subsampled GoM domain
% wrt the full TSIS domain
				ln1s = min(min(LONH));
				lt1s = min(min(LATH));
				ln2s = max(max(LONH));
				lt2s = max(max(LATH));
				D=sqrt((LON-ln1s).^2+(LAT-lt1s).^2);
				[joff,ioff]=find(D==0);

			end
%
% Subtract spatial mean:
      dmm=FF;
%      dmm(INH==0)=nan;
%      tMh=nanmean(nanmean(dmm));
      dmm(1:121,312:end)=nan; % Caribbean 
      dmm(:,504:end)=nan;
      tMh = nanmean(nanmean(dmm));
      thm = dmm-tMh;
%      thm(INH==0)=nan;


      f_pltT = 0;
      if f_pltT==1
        figure(11); clf; 
        hold on
        pcolor(LONN,LATN,tnm); shading flat;
        contour(LONN,LATN,tnm,[2.5 2.5],'k-');
        axis('equal');
        colorbar
        caxis([-4 4]);
        stl = sprintf('NEMO, dT z=%i <T>=%5.2f %i/%i/%i',abs(Z0),tMn,dm,mo,yr);
        title(stl);
        bottom_text(btx,'pwd',1);

        figure(12); clf;
        hold on;
        pcolor(LONH,LATH,thm); shading flat;
        contour(LONH,LATH,thm,[2.5 2.5],'k-');
        axis('equal');

        colorbar
        caxis([-4 4]);
        stl = sprintf('HYCOM %s, dT z=%i <T>=%5.2f %i/%i/%i',fcstname,abs(Z0),tMh,dm,mo,yr);
        title(stl);
        bottom_text(btx,'pwd',1);
%keyboard
      end       
 

      RMSE=[];
      cnn=0; 
      for inx=1:Nocn
        if mod(inx,10000)==0,
          fprintf(' Calculating RMSE: %4.1f%% done ...\n',inx/Nocn*100);
        end

        ih0=Iocn(inx);  % HYCOM index
        [j0,i0] = ind2sub(size(HH),ih0);
        ln0=LON(j0,i0);
        lt0=LAT(j0,i0);
        if ln0<ln1s | ln0>ln2s | ...
           lt0<lt1s | lt0>lt2s
          ih0=-999;  % outside subsampled GoM deep region
        end        

        II=find(IJhycom==ih0);
        if isempty(II); 
          RMSE(inx)=nan;
          continue; 
        end; % pnt outside saved indices

        Inemo=INDX.I_NEMO(II,:);
        Jnemo=INDX.J_NEMO(II,:);

        f_chck=0;
        if f_chck==1
          figure(10); clf;
          hold on;
          plot(LON(ih0),LAT(ih0),'o');

          for itt=1:4
            if itt>1
            plot(LONN(Jnemo(itt),Inemo(itt)),LATN(Jnemo(itt),Inemo(itt)),'k*');
            else
            plot(LONN(Jnemo(itt),Inemo(itt)),LATN(Jnemo(itt),Inemo(itt)),'r*');
            end
          end
        end

%
% Take the smallest RMSD
% out of the surrounding closest NEMO points:
        Err = [];
        ErrP= [];
        jhm = j0-joff+1;
        ihm = i0-ioff+1;
        if abs(HH(j0,i0)-HHx(jhm,ihm))>1e-6
          error('Check indexing for full TSIS and GoM subsampled domains ');
        end

        eTh = thm(jhm,ihm);
        for itt=1:4
          in=Inemo(itt);
          jn=Jnemo(itt);
          Err(itt)=(eTh-tnm(jn,in)).^2;
          ErrP(itt)=(tnmP(jn,in)-tnm(jn,in)).^2; % persistence
        end

        cnn=cnn+1;
        RMSE(inx)=min(Err);
        RMSEP(inx)=min(ErrP);
      end;  % ocean points
      rms=sqrt(nansum(RMSE)/cnn);
      rmsp=sqrt(nansum(RMSEP)/cnn);

% Spatial map of RMSE;
      f_chck=0;
      if f_chck==1
        ERR=zeros(m,n)*nan;
        ERR(Iocn)=sqrt(RMSE);
        figure(11); clf;
        pcolor(ERR); shading flat;
        set(gca,'xlim',[1 600],...
                'ylim',[380 850]);
        colorbar

        keyboard
      end

      RMSE=RMSE(:);
      RMSERR.ERR_squared(:,ii)=RMSE;
      RMSERR.ERRprst_squared(:,ii)=RMSEP;  % persistence squared err
      RMSERR.MeanGoM_T_hycom(ii,1) = tMh; 
      RMSERR.MeanGoM_T_nemo(ii,1)  = tMn;

      fprintf('1 record %6.3f min, RMSE=%6.4g, PerstRMSE=%6.4g, MnT=%4.1f MnT_nemo=%4.1f \n',...
              toc/60,rms,rmsp,tMh,tMn);

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

