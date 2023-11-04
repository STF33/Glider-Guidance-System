    function cmpF = smooth_colormap(cmp,nav,varargin);
%
%    function cmpF = smooth_colormap(cmp,nav);
% To make colormap look smoother without
% distinct transitions between different 
% colors, use smoother 
% that is a simple averaging
% Input: cmp - Nx3 colormap RGB
%        nav - averaging window, # of itnervals
% can specify weight for smoothing: A0=max weight
% if A0=1 or not specified - simple mean is used
% the higher A0 the more sharp colorbar is
% Note if A0 = 1 end colors are not preserved -
% they get blended with neighbor colors
%
% Dmitry Dukhovskoy, COAPS FSU
% 2014
% 2017 - added options via varargin
%
if nav<2
  fprintf('smooth_colormap: nav has to be>1, nav=%i\n',nav);
  fprintf('smooth_colormap: no smoothing is done ...\n');
  cmpF=cmp;
  return
end


nint=size(cmp,1);
if isempty(nav);
  nav=round(0.15*nint);
end
%keyboard

A0=1;
if nargin==1
  A0=varargin{1}; % weighted averaging
end;

nav2=floor(nav/2);
dA0=A0/(nav2+1);
wlft=dA0*[1:nav2]';
wrght=flipud(wlft);
wght=[wlft;A0;wrght];
lw=length(wght);

for ik=1:nint
  i1=ik-floor(nav/2);
  i2=ik+floor(nav/2);
  i1=max([i1,1]);
  i2=min([i2,nint]);
  dL=ik-i1;
  is=nav2-dL+1;
  dR=i2-ik;
  ie=nav2+1+dR;
  W=wght(is:ie);  
  for cl=1:3
    aa=cmp(i1:i2,cl);
    waa=sum(aa.*W)./sum(W);
    cmpF(ik,cl)=waa;
  end
  
%  amm=sum(cmp(i1:i2,:));
%  cmpF(ik,:)=amm;
end

return