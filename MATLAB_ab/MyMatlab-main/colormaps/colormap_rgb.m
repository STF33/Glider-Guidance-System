     function  cmp = colormap_rgb(nint);
%
% function cmp = colormap_rgb(nint);
%
% smooth transition from dark blue to green to red - white
% the color table is created by 
% creating a 3D matrix: Red changes from 0 to 1 by dR
% for each dR, G and B varies from 0 to 1
%
% nint - optional, # of intervals
% min nint = 16
% 

if nargin == 0, nint=64; end;

nR


