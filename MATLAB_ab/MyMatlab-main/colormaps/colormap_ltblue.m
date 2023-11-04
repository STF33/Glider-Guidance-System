  function cmp = colormap_ltblue(nint);
%
% creates light blue to dark skye blue
% nint - optional, # of intervals
% 
if nargin == 0, nint=10; end;


c1=zeros(nint,1);
c2=zeros(nint,1);
c3=zeros(nint,1);

p=255;
c1s=244;
c1e=0;
c2s=255;
c2e=110;
c3s=255;
c3e=255;
dc1=(c1e-c1s)/(nint-1);
dc2=(c2e-c2s)/(nint-1);
dc3=(c3e-c3s)/(nint-1);
%c2e=c2s+dc2*nint;

%c1=(c1s:dc1:c1e)'/p;
%c2=(c2s:dc2:c2e)'/p;
c1=((0:nint-1)*dc1+c1s)'/p;
c2=((0:nint-1)*dc2+c2s)'/p;
c3=((0:nint-1)*dc3+c3s)'/p;

c3(c3<0)=0;
c2(c2<0)=0;
c1(c1<0)=0;

cmp=[c1,c2,c3];

return
