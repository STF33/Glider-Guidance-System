   function [UR,VR] = rotate2XYpang(PANG, uu, vv);
%
% In netcdf HYCOM
% Output U and V components are true east/north components
% need to orient them in Cartesian coordinates, XY
% Rotate U relative to XY-grid of the figure
% Since this is curvilinear grid, at each grid point
% true North and East direction are different
%
% Input: PANG - "pang" anlge at grid point (radians!), 
%        uu, vv - true East/North components
% to get pang - use read_pang.m

% See Alan's email:
% The native model velocity variables are x-wards and y-wards (your #2), 
% but as part of the interpolation to fixed vertical levels we rotate them 
% to eastwards and northwards (your #1).  
% These are the same south of 47N, 
% where the grid is recti-linear, 
% but not north of 47N, where the grid is curvi-linear 
% in the Arctic bi-polar patch part of the tri-pole grid.
%
% If you need velocities North of 47N you have two options:
%
%  a) Use the GLBu0.08 grid fields, 
%  which are rectilinear everywhere but do not go north of 80N.
% 
% b) Rotate the velocities back to x-wards y-wards.  
%The array pang in regional.grid can be used to do this.  
%The eastwards,northwards to x-wards,y-wards code is in ALL/force/src/wi.f:
%
%              COSPANG  = COS(PANG(I,J))
%              SINPANG  = SIN(PANG(I,J))
%              TXM(I,J) = COSPANG*TXMIJ + SINPANG*TYMIJ
%              TYM(I,J) = COSPANG*TYMIJ - SINPANG*TXMIJ
%

%Unr=[uu;vv];
%Rot=[cos(lmb), sin(lmb); -sin(lmb),cos(lmb)];
%Ur=Rot*Unr;

cospang=cos(PANG);
sinpang=sin(PANG);
UR = cospang.*uu + sinpang.*vv;
VR = cospang.*vv - sinpang.*uu;

return
