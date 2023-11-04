% Plot interpolated T onto z level
% interpolation performed in: interp_hycom2z.m
%
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
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_TSIS/interp2grid
startup;

close all
clear

Z0 = -200;  % depth
irec = 4; 

EXON = zeros(9,1);
EXON(9) = 1; % select expt to plot
%ixx = 3;  % hindcast/free run # expt



pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

% HYCOM-TSIS hindcast experiments:
fhnd = 'hycom_tsis_expts.mat';
load(fhnd);
Nruns = length(EXPT);

for ii=1:Nruns
  if EXON(ii)==0
    fprintf('%i : OFF    %s \n',ii,EXPT(ii).Name);
  else
    fprintf('%i : ON ---> %s \n',ii,EXPT(ii).Name);
  end
end


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

%
% HYCOM:
rg=9806;  % convert pressure to depth, m
huge=1e20;

%Read HYCOM topography:
%ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
%HH  = -1*(nc_varget(ftopo,'mdepth'));
%LAT = nc_varget(ftopo,'mplat');
%LON = nc_varget(ftopo,'mplon');
%[mh,nh]=size(HH);
%m=mh;
%n=nh;
%HH(isnan(HH))=100;


% GoM region HYCOM:
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

Iexpt = find(EXON==1);
for jj=1:length(Iexpt);
  ixx = Iexpt(jj);
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;
%  fmatout = sprintf('%shycom_LCcontour_%2.2i.mat',pthmat,ixx);
  fmatout = sprintf('%shycom_t2Z%4.4i_hindcast%2.2i.mat',pthmat,abs(Z0),ixx);

  fprintf('Loading %s\n',fmatout);
  load(fmatout);


  TM = TZH.TM;
  dnmb = TM(irec);
  fprintf('Plotting %s\n',datestr(dnmb));

  FF = squeeze(TZH.Tz(irec,:,:));

  if ~exist('INH','var');
    HHx = TZH.HH;
    [mm,nn]=size(HHx); 
    [XM,YM]=meshgrid([1:nn],[1:mm]);
    INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
    clear XM YM;
  end

  % Subtract spatial mean T
  dmm=FF;
  dmm(INH==0)=nan;
  tM=nanmean(nanmean(dmm));
  tnm = FF-tM;
  tnm(INH==0)=nan;


  figure(1); clf;
%  pcolor(FF); shading flat;
%  caxis([16 26]);
  pcolor(tnm); shading flat;  % T demeaned
  colorbar 

end

