    function CLR = mix_2colors(cl1,cl2,N);
%  function CLR = max_2colors(cl1,cl2,N);
%
% Gradually mixes 2 colors over N intervals 
% Returns RGB colorcode (N,1:3)
%  
% Dmitry Dukhovskoy, COAPS FSU, 2/15/2013
%
if nargin == 2, nint=10; end;

nint=N;

if nint<2,
  fprintf('mix_2colors: number of intervals is small %i\n',nint);

  CLR=0.5*(cl1+cl2);
  CLR(CLR>1)=1;
  CLR(CLR<1e-3)=0;
  
  return
end

% Set Gaussian for :
r0=cl1(1);
if r0==0,
  r0=1e-6;
end
rE=cl2(1);
rE=max([rE,1e-6]);

g0=cl1(2);
if g0==0,
  g0=1e-6;
end
gE=cl2(2);
gE=max([gE,1e-6]);

b0=cl1(3);
if b0==0
  b0=1e-6;
end;
bE=cl2(3);
bE=max([bE,1e-6]);


iR=1;
iG=1;
iB=1;
if r0>rE
  dmm=rE;
  rE=r0;
  r0=dmm;
  iR=-1;
end

if g0>gE
  dmm=gE;
  gE=g0;
  g0=dmm;
  iG=-1;
end

if b0>bE
  dmm=bE;
  bE=b0;
  b0=dmm;
  iB=-1;
end

rW=sqrt(-(1-nint)^2/(log(1e-6)));


for k=1:nint
  G=exp(-(k-nint)^2/(2*rW^2));
  c1(k,1)=G;
  c2(k,1)=G;
  c3(k,1)=G;
end;
%c1(1)=r0;
%c2(1)=g0;
%c3(1)=b0;
%keyboard

dc=abs(rE-r0);
c1=c1*dc+r0;
dc=abs(gE-g0);
c2=c2*dc+g0;
dc=abs(bE-b0);
c3=c3*dc+b0;

if iR<0;
  c1=flipud(c1);
end
if iG<0;
  c2=flipud(c2);
end
if iB<0;
  c3=flipud(c3);
end


CLR=[c1,c2,c3];
CLR(CLR>1)=1;
CLR(CLR<1e-3)=0;