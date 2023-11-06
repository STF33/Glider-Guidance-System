   function Hdrf = hausdorff_distance(P,Q);
%
% Compute Hausdorff distance between two polygonal curves
% to measure similarity between two curves 2-dim, 
%  P and Q
% Coordinates are geographic coordinates, (lon,lat)
% Arrays are structured as: rows locations, columns - dimension
% e.g., Npoints, size(Py)=[Npoints,2];
% This is a maximum distance of a set to the nearest point 
% in the other set
% 
%
% Normalize:
%Py=Py./max(Py);
%Qy=Qy./max(Qy);

N=size(P,2);
if N~=2, error('2nd dim should be 2'); end;

p=size(P,1);
q=size(Q,1);
% Min distance from points P to Q
for k=1:p
  x0=P(k,1);
  y0=P(k,2);
  dst=distance_spheric_coord(y0,x0,Q(:,2),Q(:,1));
  DPmin(k,1)=min(dst);
end

% Min distance from Q to P
for k=1:q
  x0=Q(k,1);
  y0=Q(k,2);
  dst=distance_spheric_coord(y0,x0,P(:,2),P(:,1));
  DQmin(k,1)=min(dst);
end
%keyboard
aa=logical(0);
if aa,
  mxP=max(DPmin);
  mxQ=max(DQmin);
  if mxP>mxQ,
    k1=find(DPmin==mxP);
    x0=P(k1,1);
    y0=P(k1,2);
    dd = abs(distance_spheric_coord(y0,x0,Q(:,2),Q(:,1))-mxP);
    I=find(dd==min(dd));
    x1=Q(I,1);
    y1=Q(I,2);
  else
    k1=find(DQmin==mxQ);
    x0=Q(k1,1);
    y0=Q(k1,2);
    dd = abs(distance_spheric_coord(y0,x0,P(:,2),P(:,1))-mxQ);
    I=find(dd==min(dd));
    x1=P(I,1);
    y1=P(I,2);
  end
  

  plot([x0 x1],[y0 y1],'g-','linewidth',2);
  plot(x0,y0,'k*')

end;

Hdrf=max([DQmin;DPmin]);












