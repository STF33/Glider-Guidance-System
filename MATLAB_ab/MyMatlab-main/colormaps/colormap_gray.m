  function cmp = colormap_gray(nint);
%
%
% creates white - gray - black
% nint - optional, # of intervals, shades
% first color - white
% 
%

if nargin == 0, nint=64; end;

c1=zeros(nint,1);
c2=zeros(nint,1);
c3=zeros(nint,1);

c1=(0:nint-1)'./(nint-1);

cmp=[c1,c1,c1];
cmp=flipud(cmp);
