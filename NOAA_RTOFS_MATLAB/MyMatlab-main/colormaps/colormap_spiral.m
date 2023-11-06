  function cmp = colormap_spiral(nint,varargin);
%
%  function cmp = colormap_spiral(nint);
% Possible options:
% function cmp = colormap_spiral(nint,'C0',[c1 c2 c3],...
%                         'Rr',1,'Cend',[r g b]);
% creates colormap by steping in 3D 
% RGB matrix
% nint - # of intervals
% c0   - color to start [c1,c2,c3];
%       one of the triplets has to be <0, this will
%       be the linear dimension of the spiral motion
%      it will change from 0 to 1
%       The component that needs to be changed from 0 to 1, has to be < 0
%       otherwise the first component c1 is taken as linear
%       dimension
%      Default is [-1, 0, 0]
%
% Colormap is created by moving along  spiral trajectory
%  with "radius" Rr
% i.e.: circling in 2D and moving forward along third component
% Let;s say z=c1 (c1<0) - color that is gradually changing from 0 to 1
% and x=c2, y=c3 - circularly changing for each step nint
% Note that centre of rotation is at 0.5,0.5
%
% Rr - specifies: (1) sense of rotation (<0 - clockwise, >0 - c/clockwise) and
%                (2) intensity of the colors, has to be >0 and <=1
%      in the code this is (fraction of ) radius of the circle in the XY plane,
%      if Rr = 1, one of the components will always be 1
%      and XY plane  is inscribed inside the circle (with center 0.5,0.5)
%      if Rr< about 0.7, than the circle is inside the XY plane
%      if Rr=0, then no changes in XY plane (two RGB colors are fixed), 
%   only along Z (third color)
%
% Cend - [R G B] - specify color at the end of the colormap
%                 by default it will be white
%
%  Here is example how to call the subprogram:
%    nint2=100;                    % number of intervals
%    cnt2=(c1:(c2-c1)/nint2:c2);  % construct array of intervals
%    Rr=-1;                       % colors of max intensity, clockwise rotation 
%    C0=[-1,0,0];                 % starting point, red dimension
%                                 % is Z axis
%    Cend=[0.8 0.5 0.1];         % end color, default - white
%    clr = colormap_spiral(nint2,'C0',C0,'Rr',Rr,'Cend',Cend);
%
%  D. Dukhovskoy, COAPS FSU, Dec. 2012
%
%    changes Nov. 2013 - added end color option
%
SNM={'C0';'Rr';'Cend'};
%keyboard
if nargin>0,
  for n=1:2:nargin-1
    aa=varargin{n};
    ch=strcmp(aa,SNM);
    I=find(ch==1);
    if isempty(I), 
      fprintf('ERR: Call colormap_spiral error');
      help colormap_spiral.m
      error('Unknown option:',SNM); 
    end;
    switch(I)
     case(1)
      C0=varargin{n+1};
     case(2)
      Rr=varargin{n+1};
     case(3)
      Cend=varargin{n+1};
    end
  end
end



 
cl = 'function_spiral';
d2r = pi/180;
r2d = 180/pi;

%C0 = [-1,0,1];
%nint=10;

if ~exist('Rr','var'); Rr=1; end;
if ~exist('C0','var'); C0=[-1, 0, 0]; end;

if Rr<0
  srot=-1;
else
  srot=1;
end;
Rr=abs(Rr);
if Rr>1; Rr=1; end;

A0=Rr*sqrt(0.5);  % fraction of the length from the center to the corner 1,1

I=find(C0<0);
if isempty(I),
  fprintf('%s: Default: R is taken as linear dimension \n',cl);
  I=1;
  C0(I)=-1;
end;
nI=find(C0>=0&C0<=1);
if length(nI)~=2, error('initial color C0 is not properly specified'); end;
z=0;
x=C0(nI(1))-0.5;  % reference to the center (0.5, 0.5)
y=C0(nI(2))-0.5;
% Determine linear step:
dz = 1/(nint-1);

% Determine circular step (angle)
% THere has to be at least one full rotation in XY plane
%alfE = 45;  % point (1,1) (or 0.5,0.5), end of rotation after all stepping
alf0  = atan2(y,x)*r2d;
if exist('Cend','var'); % use end color to determine vect
                        % orient. at the end (Z=1)
  xe=Cend(nI(1));
  ye=Cend(nI(2));
  alfE = atan2(ye,xe)*r2d;
else
  alfE = alf0;  % vector will point to the starting point at Z=1
end

DA = alfE-alf0;
%fprintf('Initial DA= %d \n',DA);
if (DA==0); DA=1e-9; end;
if sign(DA) ~= srot,
  DA = srot*(360-abs(DA));
end;
if abs(DA)<360, DA=sign(DA)*(abs(DA)+360); end
%fprintf('Adjusted DA= %d \n',DA);
dalf = DA/(nint-1);
XYZ=[];
XYZ(1,:)=[x+0.5,y+0.5,z];
%keyboard
for k=2:nint
  alf0=alf0+dalf;
  x=A0*cosd(alf0);
  y=A0*sind(alf0);
  z=XYZ(k-1,3)+dz;

  XYZ(k,:)=[x+0.5,y+0.5,z];
end;

%plot(XYZ(:,1),XYZ(:,2),'r.-');
% Remap back to RGB:
cmp(:,nI(1)) = XYZ(:,1);
cmp(:,nI(2)) = XYZ(:,2);
cmp(:,I) = XYZ(:,3);
%cmp(end,:)=[1,1,1];
cmp(cmp>1)=1;
cmp(cmp<0)=0;

% CHange the end color is specified:
if exist('Cend','var')
  icc=round(0.8*length(cmp));
  cl1=cmp(icc,:);
  cl2=Cend;
  ni=length(cmp)-icc;
  clrM=mix_2colors(cl1,cl2,ni);
  cmp(icc+1:end,:)=clrM;
end

return
 







  
