   function dist = distance_spheric_coord(xla1,xlo1,xla2,xlo2)
%
%   function dist = distance_spheric_coord(LAT1,LON1,LAT2,LON2)
%   units = m
%  this procedure calculates the great-circle distance between two
%  geographical locations on a spheriod given it
%  lat-lon coordinates with its appropiate trigonometric
%  signs. 
%  INPUT: xla1, xlo1 - first point coordinates (latitude, longitude)
%         xla2, xlo2 - second point
% all input coordinates are in DEGREES: latitude from 90 (N) to -90,
% longitudes: from -180 to 180 or 0 to 360,
% LAT2, LON2 can be either coordinates of 1 point or N points (array)
% in the latter case, distances from Pnt 1 (LAT1,LON1) to all pnts (LAT2,LON2)
% are calculated 
%  OUTPUT - distance (in m)
%  R of the earth is taken 6371.0 km
% 
%
%
% FSU COAPS, Dmitry Dukhovskoy
% Vincenty formula is used

if (abs(xla1)>90); error('Latitude 1 is > 90'); end;
if (abs(xla2)>90); error('Latitude 2 is > 90'); end;

R=6371.0e3;
cf=pi/180;
phi1=xla1*cf;
phi2=xla2*cf;
lmb1=xlo1*cf;
lmb2=xlo2*cf;
dphi=phi2-phi1;
dlmb=lmb2-lmb1;
%keyboard
% Central angle between 2 pnts:
dmm1=(cos(phi1).*sin(dlmb)).^2;
dmm2=(cos(phi2).*sin(phi1)-sin(phi2).*cos(phi1).*cos(dlmb)).^2;
dmm3=abs(sin(phi2).*sin(phi1)+cos(phi2).*cos(phi1).*cos(dlmb));

dsgm = atan(sqrt((dmm1+dmm2)./dmm3));

% The great-circle distance:
dist1 = R*dsgm;

% Another formula:
dmm = sin(phi1).*sin(phi2)+cos(phi1).*cos(phi2).*cos(lmb2-lmb1);
dmm(dmm>1)=1;
dist2 = acos(dmm).*R;

dist=0.5*(dist1+dist2);

%keyboard

return




