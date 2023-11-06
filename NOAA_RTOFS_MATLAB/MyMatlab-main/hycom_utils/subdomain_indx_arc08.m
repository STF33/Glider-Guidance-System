    function IND = smaller_domain_indices(regn);
% Define indices for subsampling
% Greenland region
% indices are for the whole ARCc0.08 domain
% For netCDF index 1 = ind1-1, etc.

switch lower(regn)
  case('green');
   i1=450;
   i2=1200;
   j1=325;
   j2=1200;
 case('northatl');
  i1=350;
  i2=1500;
  j1=100;
  j2=1500;
 case('baffin');
  i1=350;
  i2=749;
  j1=200;
  j2=1199;
 case('beringstr');
  i1=230;
  i2=1130;
  j1=1275;
  j2=2185;
 case('atlarc'); %Natl and AO, no pacific
  i1=250;
  i2=1590;
  j1=100;
  j2=2300;
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