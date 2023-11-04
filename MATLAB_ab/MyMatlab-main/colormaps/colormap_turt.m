   function  cmp = colormap_blue(nint);
% function cmp = colormap_turt(nint);
% turtoise from lt - dart
% nint - optional, # of intervals

if nargin == 0, nint=10; end;

cmp=[];
p=255;

rS=0;
gS=0;
bS=127;
rE=220;
gE=240;
bE=240;

% Set log scale for blue
bW=(bE-bS)/log(nint);
%rW=log(nint-1)/(log(rE-rS))

% Set Gaussian for red  and green:
g0=gS;
if gS==0,
  g0=1e-6;
end
b0=rS;
if rS==0
  b0=1e-6;
end;

gW=sqrt(-(1-nint)^2/(log(g0/gE)));
rW=sqrt(-(1-nint)^2/(log(b0/bE)));

c1=zeros(nint,1);
c2=zeros(nint,1);
c3=zeros(nint,1);


for k=1:nint
  G=bS+round(bW*log(k));
  B=bS+0.5*round(bW*log(k));
%  R=rS+(k-1)^(1/rW);
%  G=-1+exp(gW*(k-1));
%  B=-1+exp(bW*(k-1));
%  B=gE*exp(-(k-nint)^2/(5*gW^2));
  R=rE*exp(-(k-nint)^2/(5*rW^2));
%  if G<0, G=0; end;
%  if G>1, G=1; end;
%  if R<0, R=0; end;
%  if R>1, R=1; end;
%  if B<0, B=0; end;
%  if B>1, B=1; end;
  c1(k)=R;
  c2(k)=G;
  c3(k)=B;
end;



cmp=[c1,c2,c3]./p;
cmp(cmp>1)=1;
cmp=flipud(cmp);

return
