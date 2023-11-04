% Extract and save LC contours from HYCOM-TSIS and
% nemo simulations
%
% Compare LC contours from the HYCOM_TSIS hindcast
% assimilating NEMO 1/100 free running fields
% and NEMO LC contours
%
% Plot HYCOM and NEMO on 1 figure
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

%f_mat = 1; % save mat
f_mat = 2; % save mat; =2 - load saved and finish missing dates


% Hindcasts:
pthd1  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
pthd2  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
pthd12 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_obs/';  % all satellites included
pthd3  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/2011_GLfreerun/'; % free run
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

fmatout = sprintf('%sLC_coord_osse_hycom_nemo.mat',pthmat);

btx = 'extr_lc_hycom_nemo.m';

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

% GoM region HYCOM:
GOM=[366   489
   476   531
   583   560
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

LCXY         = struct; 
LCXY(1).Info = 'LC fronts, NEMO 1/100 free run';
LCXY(1).Name = 'NEMO';
LCXY(2).Info = 'LC fronts, HYCOM-TSIS OSSE hindcast';
LCXY(2).Name = 'HYCOM-TSIS ugos';
LCXY(3).Info = 'LC fronts, HYCOM-TSIS OSSE Extended hindcast';
LCXY(3).Name = 'HYCOM-TSIS ugosEXT';
LCXY(4).Info = 'LC fronts, HYCOM free run';
LCXY(4).Name = 'HYCOM freerun';

if f_mat==2
  fprintf('Loading %s\n',fmatout);
  load(fmatout);

  TM1 = LCXY(1).TM;
  keyboard
end

cntr=0;
fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';
for ii=1:nrc
  tic;
%
% NEMO data:
%ncdump -h https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc
%fin='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc';
%ncid = netcdf.open(fin);
%vid  = netcdf.inqVarID(ncid,'ssh');
%ssh  = netcdf.getVar(ncid,vid,'single');
%[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,vid)
%varname = netcdf.inqVar(ncid,vid)

  yr   = YPLT(ii,1);
  mo   = YPLT(ii,2);
  dm   = YPLT(ii,3);
  dyr  = YPLT(ii,4);
  dnmb = YPLT(ii,5);
  iday = dyr;

  dnmb1=datenum(yr,mo,1);
  dnmb2=dnmb1+32;
  v2=datevec(dnmb2);
  dnmb2=datenum(v2(1),v2(2),1);
  d2=dnmb2-datenum(yr,mo,1);


  if f_mat==2
    I=find(TM1==dnmb);
    if ~isempty(I), 
      fprintf('Already extracted, skipping %s\n',datestr(dnmb));
      cntr=cntr+1;
      continue; 
    end;
  end

  fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
  fprintf('Reading NEMO: %s\n',fnemo);

  fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);

  if ~exist('LONN','var')
    fmesh=sprintf('%smesh_mask.nc',fpwd);
    dmm = ncread(fmesh,'nav_lon');
    LONN = dmm';
    dmm = squeeze(ncread(fmesh,'nav_lat'));
    LATN = dmm';

    [mm,nn] = size(LONN);

    [XM,YM]=meshgrid([1:nn],[1:mm]);
    INN = inpolygon(XM,YM,GOMN(:,1),GOMN(:,2));
    clear XM YM

  end


  enm = squeeze(ncread(fin,'ssh',[1 1 dm],[Inf Inf 1]));
  enm = enm';
  I=find(enm==0);
  enm(I)=nan;

% Subtract spatial mean ssh
  dmm=enm;
  dmm(INN==0)=nan;
  sshM=nanmean(nanmean(dmm));
  enm = enm-sshM;

  Bisol=0.17;
  dmm=enm;
  dmm(INN==0)=nan;
  LCN = identify_LC(LONN,LATN,dmm,Bisol);

% Read in HYCOM ssh - hindcast UGOS
  sday=sprintf('%3.3i',iday);
  hr=0;
  fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
  finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
  fin=fina;

  ie = exist(fin,'file');

  if ~ie
    fprintf('Missing: %s\n',fin);
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
%  dmm(HH>-200)=nan;
  sshM=nanmean(nanmean(dmm));
  ssh=ssh-sshM;

%
% Derive LC contour:
% 
  dmm=ssh;
  dmm(INH==0)=nan; 
  LCH1 = identify_LC(LON,LAT,dmm,Bisol);

% Read HYCOM-TSIS hindcast extended

  fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd2,yr,iday);
  finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd2,yr,iday);
  fin=fina;

  ie = exist(fin,'file');

  if ~ie
    fprintf('Missing: %s\n',fin);
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
%  dmm(HH>-200)=nan;
  sshM=nanmean(nanmean(dmm));
  ssh=ssh-sshM;

%
% Derive LC contour:
% 
  dmm=ssh;
  dmm(INH==0)=nan; 
  LCH2 = identify_LC(LON,LAT,dmm,Bisol);

% -----------
% Free run - shorter than hindcasts 
  fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd3,yr,iday);
  finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd3,yr,iday);
  fin=fina;
  
  ie = exist(fin,'file');
  
  if ie
				fprintf('Reading %s\n',fina);
				fld = 'srfhgt';
				[F,nn,mm,ll] = read_hycom(fina,finb,fld);
				F(F>huge)=nan;
				ssh=squeeze(F)./(1e-3*rg);  % ssh m
		%
		% Subtract anomaly:
				dmm=ssh;
				dmm(INH==0)=nan;
		%  dmm(HH>-200)=nan;
				sshM=nanmean(nanmean(dmm));
				ssh=ssh-sshM;

		%
		% Derive LC contour:
		% 
				dmm=ssh;
				dmm(INH==0)=nan;
				LCH3 = identify_LC(LON,LAT,dmm,Bisol);
  else
    fprintf('Free run is missing %s\n',datestr(dnmb));
    LCH3=[];
  end


  cntr=cntr+1;

% Nemo
  LCXY(1).TM(cntr)    = dnmb;
  LCXY(1).XY(cntr).X  = LCN.xx;
  LCXY(1).XY(cntr).Y  = LCN.yy;
% HYCOM-TSIS
  LCXY(2).TM(cntr)    = dnmb;
  LCXY(2).XY(cntr).X  = LCH1.xx;
  LCXY(2).XY(cntr).Y  = LCH1.yy;
% HYCOM-TSIS - extended
  LCXY(3).TM(cntr)    = dnmb;
  LCXY(3).XY(cntr).X  = LCH2.xx;
  LCXY(3).XY(cntr).Y  = LCH2.yy;
% HYCOM free run
  if ~isempty(LCH3)
    LCXY(4).TM(cntr)    = dnmb;
    LCXY(4).XY(cntr).X  = LCH3.xx;
    LCXY(4).XY(cntr).Y  = LCH3.yy;
  end

  fprintf('Processed 1 rec, %6.4f min\n\n',toc/60);


  if mod(ii,30)==0 & f_mat>0
    fprintf('Saving %s\n',fmatout);
    save(fmatout,'LCXY');
  end

end

fprintf('Finished, Saving %s\n',fmatout);
save(fmatout,'LCXY');


