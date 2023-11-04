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
% Sept 2021: a bug is fixed - error distance for points in different h/spheres
% exactly at opposite sides of the Globe, i.e. dlt lat = 90 and
%   dlt long = 180 dgr
% Vincenty formula on a sphere is used

if (abs(xla1)>90); error('Latitude 1 is > 90'); end;
if (abs(xla2)>90); error('Latitude 2 is > 90'); end;


%
% This formula does not work for anitpod points on the sphere!!!
% Need to fix this - better use updated version with Lambert formula
R=6371.0e3;
cf=pi/180;
phi1=xla1*cf;
phi2=xla2*cf;
lmb1=xlo1*cf;
lmb2=xlo2*cf;
dphi=abs(phi2-phi1);
dlmb=abs(lmb2-lmb1);

% Central angle between 2 pnts:
dmm1=(cos(phi1).*sin(dlmb)).^2;
dmm2=(cos(phi2).*sin(phi1)-sin(phi2).*cos(phi1).*cos(dlmb)).^2;
dmm3=sin(phi2).*sin(phi1)+cos(phi2).*cos(phi1).*cos(dlmb);

dsgm = atan(sqrt(dmm1+dmm2)./dmm3);

% The great-circle distance:
dist1 = R*dsgm;

% Another formula:
dmm = sin(phi1).*sin(phi2)+cos(phi1).*cos(phi2).*cos(lmb2-lmb1);
dmm(dmm>1)=1;
dsgm2 = acos(dmm); % cental angle
dist2 = dsgm2.*R;

%keyboard


dd=abs(1-dist1./dist2);
if max(max(dd))>1e-3,
%  I=find(dd==max(max(dd)));
%  error(' Spheric distance mismatch: %8.6g %8.6g',dist1(I),dist2(I));
  dist = dist_lmbrt;
else
  dist=0.5*(dist1+dist2);
end

%keyboard

return




