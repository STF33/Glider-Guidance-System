% Extract NEMO T field at depth z0
% demean if dmean>0
% from NEMO 1/100 simulation 
% Make fields similar to HYCOM - only 
% interior GoM discard Caribean Sea 
%
function tnm = sub_getTz0_nemo(dnmb,iz0);

fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';

DV = datevec(dnmb);
yr = DV(1);
mo = DV(2);
dm = DV(3);
dnmb1 = datenum(yr,mo,1);
dnmb2 = dnmb1+32;
v2 = datevec(dnmb2);
dnmb2 = datenum(v2(1),v2(2),1);
d2 = dnmb2-datenum(yr,mo,1);

tnm = sub_get_NEMO_TS(dnmb,'toce',iz0);
tnm(1:908,998:end) = nan;
tnm(908:end,1732:end) = nan;

return
