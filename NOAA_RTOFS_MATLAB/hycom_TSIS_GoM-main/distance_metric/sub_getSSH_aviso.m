%   Get SSH from AVISO
% Demeaned or not

adt = squeeze(nc_varget(flaviso,'adt'));
adt = adt(j1A:j2A,i1A:i2A);

