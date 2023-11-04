     function  cmp = colormap_hot(nint);
%
% function cmp = colormap_hot(nint);
%
% creates dark red  to light red to white c/bar
%
% nint - optional, # of intervals
% 

if nargin == 0, nint=64; end;

w3=0.375;
w2=0.375;
w1=1-(w2+w3);

n3=round(nint*w3);
n2=round(nint*w2);
n1=nint-(n2+n3);

dc=1/n3;

c1=zeros(nint,1);
c2=zeros(nint,1);
c3=zeros(nint,1);

c1=(1:nint)'*1/n3;
c1(c1>1)=1;
i3=min(find(c1==1));
c2(i3+1:nint)=(1:nint-n2)*1/n2;
c2(c2>1)=1;
i2=min(find(c2==1));
c3(i2+1:nint)=(1:n1)*1/n1;

cmp=[c1,c2,c3];
cmp=flipud(cmp);

return
