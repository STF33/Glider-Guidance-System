  function cmp = colormap_purple(nint);
%
% creates light purple  to dark purple
% nint - optional, # of intervals
% 
if nargin == 0, nint=10; end;


c1=zeros(nint,1);
c2=zeros(nint,1);
c3=zeros(nint,1);

p=255;
c1s=p;
c1e=100;
c2s=225;
dc1=(c1e-c1s)/(nint-1);
dc3=dc1;
dc2=(0-c2s)/(0.6*(nint-1));
%c2e=c2s+dc2*nint;

c1=(c1s:dc1:c1e)'/p;
%c2=(c2s:dc2:c2e)'/p;
c2=((0:nint-1)*dc2+c2s)'/p;
c2(c2<0)=0;
c1(c1<0)=0;
c3=c1;

cmp=[c1,c2,c3];


