% read NEMO grid
% NEMO data:
%ncdump -h https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc
%fin='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc';
%ncid = netcdf.open(fin);
%vid  = netcdf.inqVarID(ncid,'ssh');
%ssh  = netcdf.getVar(ncid,vid,'single');
%[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,vid)
%varname = netcdf.inqVar(ncid,vid)

function [LONN,LATN,ZZ] = sub_get_NEMO_grid(dnmb);

fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';

DV = datevec(dnmb);
yr = DV(1);
mo = DV(2);

dnmb1=datenum(yr,mo,1);
dnmb2=dnmb1+32;
v2=datevec(dnmb2);
dnmb2=datenum(v2(1),v2(2),1);
d2=dnmb2-datenum(yr,mo,1);

fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
fprintf('Reading NEMO grid: %s\n',fnemo);

fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);

fmesh=sprintf('%smesh_mask.nc',fpwd);
dmm = ncread(fmesh,'nav_lon');
LONN = dmm';
dmm = squeeze(ncread(fmesh,'nav_lat'));
LATN = dmm';

[mm,nn] = size(LONN);

%[XM,YM]=meshgrid([1:nn],[1:mm]);
%INN = inpolygon(XM,YM,GOMN(:,1),GOMN(:,2));
%clear XM YM

ZZ = ncread(fin,'deptht');
ZZ = -ZZ;


return
