    function IND = smaller_domain_indices(regn);
% Define indices for subsampling
% Greenland region
% indices are for the whole ARCc0.08 domain
% For netCDF index 1 = ind1-1, etc.

switch lower(regn)
  case('green');
   i1=900;
   i2=2400;
   j1=650;
   j2=2400;
 case('northatl');
  i1=700;
  i2=3000;
  j1=200;
  j2=3000;
 case('baffin');
  i1=700;
  i2=1498;
  j1=400;
  j2=2398;
 case('beringstr');
  i1=460;
  i2=2260;
  j1=2550;
  j2=4370;
 case('atlarc'); %Natl and AO, no pacific
  i1=500;
  i2=3180;
  j1=200;
  j2=4600;
end
dj=j2-j1+1;
di=i2-i1+1;
IND.i1=i1;
IND.i2=i2;
IND.j1=j1;
IND.j2=j2;
IND.dj=dj;
IND.di=di;

return