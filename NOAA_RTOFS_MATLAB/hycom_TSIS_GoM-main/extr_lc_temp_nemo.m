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
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 1; % save mat; =2 - load saved and finish missing dates
T0 = 2.5;   % T contour anomaly - demeaned T
Z0 = -200;  % depth

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';



YPLT=[];
cc=0;
for iy=2011:2012
  for dd=1:365
    if iy==2011 & dd==1; continue; end;
    if iy==2012 & dd>182,
      break;
    end
    dnmb=datenum(iy,1,1)+dd-1;
    dv=datevec(dnmb);
    cc=cc+1;
    YPLT(cc,1)=iy;
    YPLT(cc,2)=dv(2);
    YPLT(cc,3)=dv(3);
    YPLT(cc,4)=dd;
    YPLT(cc,5)=dnmb;
  end
end

nrc=cc;


% GoM region, NEMO:
GOMN=[         100         365
         651         337
        1091         687
        1246         798
        1512         881
        1787         998
        1904        1292
        1710        1914
          23        1920
           8         748];


fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';

fmatout = sprintf('%sNEMO_dT%2.2i_Z%3.3i_contour.mat',pthmat,round(T0*10),abs(Z0));

LCXY = struct;
LCXY.Info = 'LC contour from NEMO T anom at Z0';
LCXY.Z0 = Z0;
LCXY.dT0 = T0;

cntr=0;
LONN=[];
LATN=[];
INN=[];
for ii=1:nrc
		tic;

  yr   = YPLT(ii,1);
  mo   = YPLT(ii,2);
  dm   = YPLT(ii,3);
  dyr  = YPLT(ii,4);
  dnmb = YPLT(ii,5);
  iday = dyr;


		DV = datevec(dnmb);
		dnmb1=datenum(yr,mo,1);
		dnmb2=dnmb1+32;
		v2=datevec(dnmb2);
		dnmb2=datenum(v2(1),v2(2),1);
		d2=dnmb2-datenum(yr,mo,1);

		fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
		fprintf('Reading NEMO: %s\n',fnemo);

		fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);

		if isempty(LONN)
				fmesh=sprintf('%smesh_mask.nc',fpwd);
				dmm = ncread(fmesh,'nav_lon');
				LONN = dmm';
				dmm = squeeze(ncread(fmesh,'nav_lat'));
				LATN = dmm';

				[mm,nn] = size(LONN);

				[XM,YM]=meshgrid([1:nn],[1:mm]);
				INN = inpolygon(XM,YM,GOMN(:,1),GOMN(:,2));
				clear XM YM

    ZZ = ncread(fin,'deptht');
    ZZ = -ZZ;

    dZ = abs(abs(ZZ)-abs(Z0));
    iz0 = find(dZ==min(dZ));
		end


		tnm = squeeze(ncread(fin,'toce',[1 1 iz0 dm],[Inf Inf 1 1]));
		tnm = tnm';
		I=find(tnm==0);
		tnm(I)=nan;

		% Subtract spatial mean ssh
		dmm=tnm;
		dmm(INN==0)=nan;
		tM=nanmean(nanmean(dmm));
		tnm = tnm-tM;

		dmm=tnm;
		dmm(INN==0)=nan;
 
 
	 Bisol=T0;
  npnt0=200;
  [XX,YY] = get_contour(LONN,LATN,dmm,Bisol,npnt0); 

% Nemo
		cntr=cntr+1;
		LCXY.TM(cntr)    = dnmb;
		LCXY.XY(cntr).X  = XX;
		LCXY.XY(cntr).Y  = YY;

		fprintf('Processed 1 rec, %6.4f min\n',toc/60);

		if mod(ii,30)==0 & f_mat>0
				fprintf('Saving %s\n',fmatout);
				save(fmatout,'LCXY');
		end

end

fprintf('Saving NEMO %s\n',fmatout);
save(fmatout,'LCXY');



return
