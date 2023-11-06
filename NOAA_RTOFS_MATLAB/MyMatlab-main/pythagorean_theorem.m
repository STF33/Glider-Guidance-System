   function [distlat, distlong, dist] = distance_spheric_coord(R,phi1,lmb1,...
	 phi2,lmb2)

% NOTE: THIS METHOD IS ONLY FOR SPHERICAL TRIANLGES, i.e.
% triangles formed by great circles (centered in the center of the sphere)
% it is not accurate if two points are far away from equator !!!

%   function [distx, disty, dist] = distance_spheric_coord(R,phi1,lmb1,...
%	 phi2,lmb2)
% calculates shortest distances between
% two points in spherical polar coordiantes (length
% of the shortest arc between 2 points on a circle)
% Spherical angle for computing total distance is calculated
% using spherical Pythagorean theorem
%
% Input: R - radius of the Earth
% phi - latitude (anlge between equator plane and radius vector)
% lmb - longitude
% All input angles are in degrees !!!
% 
% output: distlat -  zonal distance calculated at latitude phi1
%         distlon - meridional distance
%         dist - length of the arc connecting two points


  alf=abs(phi2-phi1)*pi/180;
	bet=abs(lmb2-lmb1)*pi/180;
	
	gam=acos(cos(alf)*cos(bet));

  rr=abs(R*cos(phi1*pi/180));
  distlat=rr*bet;
	distlong=R*alf;
	dist=rr*gam;
	
	
