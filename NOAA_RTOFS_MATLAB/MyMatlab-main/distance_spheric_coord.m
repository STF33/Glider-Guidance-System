   function dist_lmbrt = distance_spheric_coord(xla1,xlo1,xla2,xlo2)
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
%
% FSU COAPS, Dmitry Dukhovskoy
% Sept 2021: Lambert formula for distance on an ellipse

if (abs(xla1)>90); error('Latitude 1 is > 90'); end;
if (abs(xla2)>90); error('Latitude 2 is > 90'); end;

%
% This formula does not work for anitpod points on the sphere!!!
R=6371.0e3;
cf=pi/180;
phi1=xla1*cf;
phi2=xla2*cf;
lmb1=xlo1*cf;
lmb2=xlo2*cf;
dphi=abs(phi2-phi1);
dlmb=abs(lmb2-lmb1);

I=[];
[a1,a2] = size(dphi);
if a1>1 | a2>1
	I=find(dphi==0 & dlmb==0);
	dphi(I)=1.e-30;
	dlmb(I)=1.e-30;
else
	if dphi==0 & dlmb==0
		dist_lmbrt = 0;
		return
	end
end

%
% The most accurate is Lambert's formula:
% Lambert's formula is the method used to calculate the shortest distance
% along the surface of an ellipsoid. When used to approximate the Earth 
% and calculate the distance on the Earth surface, it has an accuracy on 
% the order of 10 meters over thousands of kilometers, which is more precise 
% than the haversine formula.
Req = 6378.e3;
Rpl = 6357.e3;
Eflat = (Req-Rpl)/Req; % flatenning of the Earth

% Haversine formula to calculate central angle:
aa1 = (sin(dphi/2.)).^2;
aa2 = cos(phi1).*cos(phi2).*(sin(dlmb/2.)).^2;
dsgm_hv = 2*asin(sqrt(aa1+aa2));

%
% Reduced latitudes - due to flattening:
beta1 = atan((1-Eflat)*tan(phi1));
beta2 = atan((1-Eflat)*tan(phi2));
PP = 0.5*(beta1+beta2);
QQ = 0.5*(beta2-beta1);
X = (dsgm_hv-sin(dsgm_hv)).*( (sin(PP)).^2.*(cos(QQ)).^2 )./( (cos(dsgm_hv/2).^2) );
Y = (dsgm_hv+sin(dsgm_hv)).*( (cos(PP)).^2.*(sin(QQ)).^2 )./( (sin(dsgm_hv/2).^2) );


dist_lmbrt = Req*(dsgm_hv-Eflat/2*(X+Y));
if ~isempty(I)
	dist_lmbrt(I) = 0.;
end


return




