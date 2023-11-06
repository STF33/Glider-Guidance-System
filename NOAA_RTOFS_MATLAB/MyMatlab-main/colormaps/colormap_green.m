  function cmp = colormap_green(nint);
% creates light green  to dark green
% nint - optional, # of intervals
% 
if nargin == 0, nint=10; end;


c1=zeros(nint,1);
c2=zeros(nint,1);
c3=zeros(nint,1);

p=255;
c1s=240;
c1e=0;
dc1=(c1e-c1s)/(0.5*(nint-1));
%c1e=c1s+dc1*(nint-1);


c2s=p;
c2e=30;
dc2=(c2e-c2s)/(0.5*(nint-1));

%c1=(c1s:dc1:c1e)'/p;
c1=((0:nint-1)*dc1+c1s)'/p;
%c2=(c2s:dc2:c2e)'/p;
c2=((0:nint-1)*dc2-(nint-1)*dc2+c2e)'/p;
c2(c2<0)=0;
c2(c2>0.99)=0.99;
c1(c1<0)=0;
c1(c1>0.99)=0.99;
c3=c1;

cmp=[c1,c2,c3];






