% Predictability forecasts
% For LC contour finding: interpolate T or S fields to 
% Z depth Z0
%
% Extract LC and LCEs for calculating MHD
% Using T fields (demeaned) at depth Z0, contour T anomaly
%
% NEMO data
%
% NEMO data:
%ncdump -h https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc
%fin='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc';
%ncid = netcdf.open(fin);
%vid  = netcdf.inqVarID(ncid,'ssh');
%ssh  = netcdf.getVar(ncid,vid,'single');
%[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,vid)
%varname = netcdf.inqVar(ncid,vid)
%
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_TSIS/interp2grid
startup;

close all
clear

f_mat = 1; % save mat; =2 - load saved and add  missing dates
Z0 = -200;  % depth
iFcst = 16;  % #10, ..., 16

if f_mat>2; fprintf('wrong f_mat flag %i\n',f_mat); return; end;

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';


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

%
% Indices for subsampled GoM domain 
[XX,YY,HHs,INH,Iocn,xt1,xt2,yt1,yt2,its1,its2,jts1,jts2] = ...
   sub_subsample_HYCOM2GoM(HH,LON,LAT,Z0);

Hnas=HHs;
[msb,nsb]=size(HHs);

dz=50;
ZMnas = [Z0+dz:-dz:Z0-dz];
iz0 = find(ZMnas==Z0);

knas=length(ZMnas);
%[mnas,nnas]=size(nlon);
mnas = m;
nnas = n;

% Create Land/bottom mask indices:
for kk=1:knas
		I=find(Hnas>=ZMnas(kk));
		ILAND(kk).I=I;
end

% 
% 1 forecast group at a time
Nhind = iFcst;
hnd_name = FCST.Hindcast_Name;
pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
irun1 = FCST.run1;
irun2 = FCST.run2;
ntime = FCST.ntime;
for itime=1:ntime
  for irun=irun1:irun2
    clear TZH
				TZH.Info = 'Interpolated T to fixed z0';
				TZH.Z0   = Z0;
				TZH.HH   = HHs;
				TZH.LON  = XX;
				TZH.LAT  = YY;
				TM = [];

				RUN = FCST.TIME0(itime).RUN(irun);
				pthd1 = RUN.pthbin;
				TM    = RUN.TM;
				YDAY  = RUN.jday;
				nrc   = length(TM);
				DV    = datevec(TM);

    nmexp = sprintf('fcst%2.2i-%2.2i%2.2i',Nhind,itime,irun);

    TZH.Name = nmexp;
    TZH.Pthdata = pthd1;

    fmatout = sprintf('%shycom_t2Z%4.4i_fcstPrdct%2.2i-%2.2i%2.2i.mat',...
                    pthmat,abs(Z0),Nhind,itime,irun);
				fprintf('\n\n %s %s\n',nmexp,fmatout);
    fprintf('Forecast Name: %s\n',hnd_name);
				fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',Nhind,itime,irun);
				fprintf(' Input data: %s\n',pthd1);


    ilast=0;
    TM_saved=[];
				if f_mat==2
      if exist(fmatout,'file')
  						fprintf('\n\n !!! Continue from last saved record, loading %s\n',fmatout);
		  				load(fmatout);
  						TM_saved = TZH.TM;
        f_save=0;  % do not save if whole run is already extracted
       
%
% Find last record in mat file
        ilast=length(TM_saved);
% Check that the dates are correct:
        if ilast>0 
          dT = TM_saved-TM(1:ilast); 
          if max(dT)>0
            error('Saved mat file - wrong dates');
          end
        end
      else
        fprintf('%s not found to continue, start from day 1 %s\n',fmatout,datestr(TM(1)));
        ilast=0;
      end
				end
    if f_mat==0, f_save=0; end;

				irec = 0;
				for ii=1:nrc
						yr  = DV(ii,1);
						mo  = DV(ii,2);
						dm  = DV(ii,3);
						dnmb= TM(ii); 
						iday= YDAY(ii);

      if isempty(dnmb); keyboard; end;
      if ii<=ilast
        fprintf('irec=%i, Date %s found, skipping\n',ii,datestr(dnmb));
        ii=ii+1;
        continue;
      else
        if f_mat>=1; f_save=1; end;
      end

%keyboard
      fprintf('Day #%i, Fcst name %s\n',ii,nmexp);
						fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
						finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
						fin=fina;

						ie = exist(fin,'file');

						if ~ie
								fprintf('  ERR ** Missing forecast: %s\n',fin);
								keyboard
								continue;
						end

		%keyboard
						tic;
						[F,n,m,l] = read_hycom(fina,finb,'thknss');
						F=F/rg;
						F(F>1e10)=0;
						dH=squeeze(F); % note this is not U,V thickness, center grid
										% there are thck-u & thck-v 
						ih=find(dH<1e-1);
						dHs=dH(:,jts1:jts2,its1:its2);

						[ZM1,ZZ1] = sub_zz_zm(fina,finb,HH,'thknss',dH);
						ZM1(isnan(ZM1))=-100;
						ZZ=ZZ1(:,jts1:jts2,its1:its2);
						ZM=ZM1(:,jts1:jts2,its1:its2);

		%keyboard
						fprintf('Reading %s\n',fina);
						fld = 'temp';
						[F,nn,mm,ll] = read_hycom(fina,finb,fld);
						F(F>huge)=nan;
						T=F;
		%    T=F(:,jts1:jts2,its1:its2);

      f_saveIFL=0;
						if ~exist('IFL','var');
								IFL=[];
								fifl = sprintf('%sIFL.mat',pthmat);
        if exist(fifl,'file')
  								load(fifl);
        else
          f_saveIFL=1;
        end
						end

		% fill nans in the surface layer  
						t1=squeeze(T(1,:,:));
%						[t1,IFL]=sub_fill_land(t1,IFL);
						t1=sub_fill_land(t1);
		%   fifl = sprintf('%sIFL.mat',pthmat);
      if f_saveIFL==1
  		    save(fifl,'IFL');
      end
						T(1,:,:)=t1;
						T0=T;
		% Subsample domain:
						T=T0(:,jts1:jts2,its1:its2);

		% Fill nans in the subsurface layers (bottom):
						T = sub_fill_3z(T);

		% Interpolate onto NAS grid:    
		% Interpolate into z levels
		% as a 2D vert sections
						[msb,nsb]=size(XX);
						Tzi=[];
						for iis=1:nsb
								if mod(iis,100)==0
										fprintf('%s interp to z, done %5.2f%% ...\n',fld,iis/nsb*100);
								end

								y1=YY(:,iis);
								ZMh=squeeze(ZM(:,:,iis));
								F=squeeze(T(:,:,iis));
								Hs = HHs(:,iis);
								Hs(Hs>=0)=-1;

								Fi=sub_interp2z_2D(Hs,F,ZMnas,ZMh);
								Tzi(:,:,iis)=Fi;

						end

		% Put land back
						for kk=1:knas
								il=ILAND(kk).I;
								Tzi(kk,il)=nan;
						end
						Tz0 = squeeze(Tzi(iz0,:,:));

		% Check:
		%    T=T0(:,jts1:jts2,its1:its2);
		%    aa = squeeze(T(13,:,:));

						irec=irec+1;
						TZH.TM(irec) = dnmb;
						TZH.Tz(irec,:,:) = Tz0;

						fprintf('min/max Tintrp= %8.2f %8.2f\n',min(min(Tz0)),max(max(Tz0)));    
						fprintf('Processed 1 rec, %6.4f min\n\n',toc/60);

						if mod(irec,15)==0 & f_save==1
								fprintf('Saving %s\n',fmatout);
								save(fmatout,'TZH');
						end

				end
%    keyboard
    if f_save==1
      fprintf('---- End ----   Saving %s\n',fmatout);
      save(fmatout,'TZH');
    else
      fprintf('No output is saved, f_mat=%i, f_save=%i\n',f_mat,f_save);
    end

  end


end

