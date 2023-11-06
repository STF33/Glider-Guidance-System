   function Hdrf = hausdorff_distance3(P,Q,varargin);
%
%  hausdorff_distance(P,Q) - calc. Hausd. dist, P,Q are 
%                            geogr. coord.
%
%  hausdorff_distance(P,Q,'geo') - calc. Hausd. dist, P,Q are 
%                            geogr. coord.
% Compute Hausdorff distance between two polygonal curves
% to measure similarity between two curves 3-dim, 
%  P and Q
% Default Coordinates are Cartesian, i.e.
% GEOGRAPHIC coordinates will be treated as (x,y)
% Options: 
%          'geogr' - lon/lat are converted into cartesian km x,y
% Euclidian distance is computed
% Arrays are structured as: rows locations, columns - dimension
% e.g., Npoints, size(Py)=[Npoints,2];
% This is a maximum distance of a set to the nearest point 
% in the other set
% 
% Option 'ipdm' - use ipdm (inter-point distance matrix) 
%                 to find Hausd. distance
%                 >6 times faster routine
%                 but may not be correct for geographic
%                 coordinates and units are meaningless
%                 
%  hausdorff_distance(P,Q,'cart','ipdm') - use cart. coord
%            calculate Hausd. using ipdm function
%
% Dmitry Dukhovskoy, 2011
% Changes: 
% 2013 - distance for both Geogr./cart. coord
%      - added ipdm option, much faster calculation
% 2014 - code is cleaned a bit
%
%keyboard
fgeogr=logical(0);
fipdm=logical(0);
cf=1e-3;         % convert m->km 
if nargin>=3
  for ik=3:nargin
    aa=varargin{ik-2};
    i1=strmatch('geo',lower(aa));
  %  i2=strmatch('cart',lower(aa));
    i2=strmatch('ipdm',lower(aa));
    if ~isempty(i1);
      fgeogr=logical(1);
    end

    if ~isempty(i2)
      fipdm=logical(1);
    end
  end
end

%fipdm=logical(0);
%if nargin==4
%  aa=varargin{2};
%  if strcmp(lower(aa),'ipdm');
%    fipdm=logical(1);
%  end
%end

% Convert to cartesian
if fipdm & fgeogr
  x0=P(1,1);
  y0=P(1,2);
  ix=find(P(:,1)<x0);
  iy=find(P(:,2)<y0);
  px=distance_spheric_coord(y0,x0,y0,P(:,1))*cf;  % m-> km
  py=distance_spheric_coord(y0,x0,P(:,2),x0)*cf;
  px(ix)=-px(ix);
  py(iy)=-py(iy);
  ix=find(Q(:,1)<x0);
  iy=find(Q(:,2)<y0);
  qx=distance_spheric_coord(y0,x0,y0,Q(:,1))*cf;  % m->km
  qy=distance_spheric_coord(y0,x0,Q(:,2),x0)*cf;
  qx(ix)=-qx(ix);
  qy(iy)=-qy(iy);
  P=[px,py];
  Q=[qx,qy];
%  keyboard
end

N=size(P,2);
if N~=3, error('2nd dim should be 3'); end;

% Check: 

p=size(P,1);
q=size(Q,1);

if fipdm
%  tic
  idst = ipdm(P,Q);
  DPmin=min(idst,[],2); % min distances P to Q
  DQmin=min(idst,[],1); % min distance Q to P
  DQmin=DQmin';  
%  toc
else
% Min distance from points P to Q
%tic
  for k=1:p
    x0=P(k,1);
    y0=P(k,2);
    z0=P(k,3);
%    if fgeogr,
%      dst=distance_spheric_coord(y0,x0,Q(:,2),Q(:,1))*cf;
%    else
%    dst=sqrt((Q(:,1)-x0).^2+(Q(:,2)-y0).^2);
    dst=sqrt((Q(:,1)-x0).^2+(Q(:,2)-y0).^2+(Q(:,3)-z0).^2);
%    end
    DPmin(k,1)=min(dst);
  end

  % Min distance from Q to P
  for k=1:q
    x0=Q(k,1);
    y0=Q(k,2);
    z0=Q(k,3);
%    if fgeogr,
%      dst=distance_spheric_coord(y0,x0,P(:,2),P(:,1))*cf;
%    else
    dst=sqrt((P(:,1)-x0).^2+(P(:,2)-y0).^2+(P(:,3)-z0).^2);
%    end
    DQmin(k,1)=min(dst);
  end

  %keyboard
  % Example - works only for geogr. coord
  aa=logical(0);
  if aa,
    mxP=max(DPmin);
    mxQ=max(DQmin);
    if mxP>mxQ,
      k1=find(DPmin==mxP);
      x0=P(k1,1);
      y0=P(k1,2);
      dd = abs(distance_spheric_coord(y0,x0,Q(:,2),Q(:,1))*cf-mxP);
      I=find(dd==min(dd));
      x1=Q(I,1);
      y1=Q(I,2);
    else
      k1=find(DQmin==mxQ);
      x0=Q(k1,1);
      y0=Q(k1,2);
      dd = abs(distance_spheric_coord(y0,x0,P(:,2),P(:,1))*cf-mxQ);
      I=find(dd==min(dd));
      x1=P(I,1);
      y1=P(I,2);
    end


    plot([x0 x1],[y0 y1],'g-','linewidth',2);
    plot(x0,y0,'k*')

  end;
%toc
end;  % if fipdm

Hdrf=max([DQmin;DPmin]);

return











