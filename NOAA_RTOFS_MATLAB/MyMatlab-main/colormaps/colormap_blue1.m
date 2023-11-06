  function cmp = colormap_blue;
% function cmp = colormap_blue(nint);
%
% creates light blue - to white c/bar
%
% nint - optional, # of intervals

disp('Not finished, use colormap_cold');
error('Not finished ...')

if nargin == 0, nint=64; end;

r1=0.75;
r2=0;

g1=0.95;
g2=0;

b1=1;
b2=0.25;

w3=(b1-b2)/(nint-1);
w2=(g1-g2)/(nint-1);
w1=(r1-r2)/(nint-1);


c3=(0:nint-1)'*w3+b2;
c3(c3>1)=1;
c2=(0:nint-1)'*w2+g2;
c2(c2>1)=1;
c1=(0:nint-1)'*w1+r1;
c1(c1>1)=1;

cmp=[c1,c2,c3];
cmp=flipud(cmp);


