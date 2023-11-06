% Get info aboutNEMO grid
%
%
function [LONN,LATN,INN] = sub_getNEMO_info;

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


yr   = 2011;
mo   = 1;
dm   = 10;
dyr  = 10;
dnmb = datenum(yr,mo,dm);
iday = dyr;


DV = datevec(dnmb);
dnmb1=datenum(yr,mo,1);
dnmb2=dnmb1+32;
v2=datevec(dnmb2);
dnmb2=datenum(v2(1),v2(2),1);
d2=dnmb2-datenum(yr,mo,1);


fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';
fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
fprintf('Reading NEMO: %s\n',fnemo);

fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);
%keyboard
		fmesh=sprintf('%smesh_mask.nc',fpwd);
		dmm = ncread(fmesh,'nav_lon');
		LONN = dmm';
		dmm = squeeze(ncread(fmesh,'nav_lat'));
		LATN = dmm';

		[mm,nn] = size(LONN);

		[XM,YM]=meshgrid([1:nn],[1:mm]);
		INN = inpolygon(XM,YM,GOMN(:,1),GOMN(:,2));



return
