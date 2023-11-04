   function  cmp = colormap_red(nint);
% function cmp = colormap_red(nint);
%
% nint - optional, # of intervals
% min nint = 16

if nargin == 0, nint=9; end;

cmp=[];
p=255;

rS=127;
gS=0;
bS=0;
rE=255;
gE=250;
bE=245;

% Set log scale for red
rW=(rE-rS)/log(nint);
%rW=log(nint-1)/(log(rE-rS))

% Set Gaussian for blue and green:
g0=gS;
if gS==0,
  g0=1e-6;
end
b0=bS;
if bS==0
  b0=1e-6;
end;

gW=sqrt(-(1-nint)^2/(log(g0/gE)));
bW=sqrt(-(1-nint)^2/(log(b0/bE)));

c1=zeros(nint,1);
c2=zeros(nint,1);
c3=zeros(nint,1);


for k=1:nint
  R=rS+round(rW*log(k));
%  R=rS+(k-1)^(1/rW);
%  G=-1+exp(gW*(k-1));
%  B=-1+exp(bW*(k-1));
  G=gE*exp(-(k-nint)^2/(5*gW^2));
  B=bE*exp(-(k-nint)^2/(5*bW^2));
  c1(k)=R;
  c2(k)=G;
  c3(k)=B;
end;



cmp=[c1,c2,c3]./p;
cmp(cmp>1)=1;
cmp=flipud(cmp);

return
