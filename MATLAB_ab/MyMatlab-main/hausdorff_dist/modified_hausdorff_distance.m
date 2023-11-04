   function  MHD  = modified_hausdorff_distance( P,Q,varargin )
%
% This function computes the modified Hausdorff distance between two sets
% P,Q entered as n x k matrices with dimension k constant (the sets can
% have different numbers of points, however). It improves over its
% predecessor by being much less sensitive to outlier points.  The function
% is defined by
%
% d_MH(P,Q) = max{1/|P|\sum_{p\in P}d(p,Q),1/|Q|\sum_{q\in Q}d(q,P)}
%
% where d(p,Q) = min_{q\in Q}d(p,q), and similarly for d(P,q).
%
% The function is not a true topological metric since it fails the triangle
% inquality requirement.  In practice, this does not appear to be a
% problem.  Indeed, the function corresponds strongly with the human
% perception of shape.
%
% Reference: "A Modified Hausdorff Distance for Object Matching," Dubuisson
% & Jain.  Proc. International Conference on Pattern Recognition,
% Jerusalem, Israel.  1994.
%
%  Options:   if P,Q are geogr. coord, its better to convert
%             them into cartesian x,y (km)
%
% add option 'geo': modified_hausdorff_distance( P,Q,'geo' )
% otherwise no conversion is done
%
% FSU COAPS
% Dmitry Dukhovskoy
% Code written with Jonathan Ubnoske
%

%keyboard
if nargin<3
  fgeogr=logical(0);
else
  aa=varargin{1};
  i1=strmatch('geo',lower(aa));
  if ~isempty(i1);
    fgeogr=logical(1);
  else
    fgeogr=logical(0);
    fprintf('\n !!! coordinate specification is not clear %s\n',aa);
    fprintf('!!! set default: not geo \n\n');
  end
end



sP = size(P);
sQ = size(Q);

% Check validity of inputs
if sP(2)~=sQ(2)
    error('Sets P and Q must be of the same dimension')
end

%keyboard

% Convert geographic -> cartesian (km)
if fgeogr
  x0=P(1,1);
  y0=P(1,2);
  ix=find(P(:,1)<x0);
  iy=find(P(:,2)<y0);
  px=distance_spheric_coord(y0,x0,y0,P(:,1))*1e-3;  % m-> km
  py=distance_spheric_coord(y0,x0,P(:,2),x0)*1e-3;
  px(ix)=-px(ix);
  py(iy)=-py(iy);
  ix=find(Q(:,1)<x0);
  iy=find(Q(:,2)<y0);
  qx=distance_spheric_coord(y0,x0,y0,Q(:,1))*1e-3;  % m->km
  qy=distance_spheric_coord(y0,x0,Q(:,2),x0)*1e-3;
  qx(ix)=-qx(ix);
  qy(iy)=-qy(iy);
  
  if sP<=2
    P=[px,py];
    Q=[qx,qy];
  else
    nE=sP(2);
    P=[px,py,P(:,3:nE)];
    Q=[qx,qy,Q(:,3:nE)];
  end
  
%  keyboard
end
%keyboard


% Vectorize.
store_dist = ipdm(P,Q);

dist_p_Q = min(store_dist,[],2);
dist_P_q = min(store_dist,[],1);

dist_P_Q = 1/sP(1)*sum(dist_p_Q);
dist_Q_P = 1/sQ(1)*sum(dist_P_q);

MHD = max(dist_P_Q,dist_Q_P);

%keyboard

f_chck=0;
if f_chck>0
  p=length(P);
  q=length(Q);
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
  dPQ=1/sP(1)*sum(DPmin);
  dQP=1/sP(1)*sum(DQmin);
  mhd2=max([dPQ,dQP]);
  
end  

return

