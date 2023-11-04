  function  [AX,AY]=runge_kutta(U1,V1,dt,...
                           U2,V2, LX, LY, DX, DY,HH,nu);
%
%  function  [AX,AY]=runge_kutta(U1,V1,dt,...
%                           U2,V2, LX, LY, DX, DY,HH,nu);
%
% Advect particle using Runge-Kutta scheme
% over 1 time step = dt
%
% Arrays: [m,n]: U1,U2,V1,V2,HH
% Input U1,V1 - u,v components of velocity fields at time level n
%       U2,V2, - u,v components at next time level n+1
%              make sure that U and V are rotated into i,j
%              so that U>0 is to the right everywhere on the grid
%              V>0 is upward everywhere
%      LX, LY - particle coordinates at time step N (np,1)
%               use indices as coordinates
%      DX, DY  - spacing in X and Y directions of the grid
%
%      HH  - topography
%      nu  - diffusion, m2/s, =0 - no diffusion
%            parameterized by Random walk
%
% Output: AX, AY - Updated coordinates of particles at time n+1
%
% U, V are X-, Y-ward components
%
% In HYCOM: if use output interpolated into fixed depths (z): 
% Rotate U relative to XY-grid of the figure
% Since this is curvilinear grid, at each grid point
% true North and East direction are different
% Find true NE direction and plot vector
%           u>0 is to the East (positive X direction)
%           v>0 is to the North (+ y direction)
% Project vector on local North and Eastward gradient vectors
%
% Use PANG angle - HYCOM output - to do the rotation 
% 
% D. Dukhovskoy, COAPS FSU, 2014
%
if size(HH)~=size(U1) |...
   size(HH)~=size(V1) |...
   size(HH)~=size(U2) |...
   size(HH)~=size(V2),
   error('Input Arrays sizes do not match')
end

[m,n]=size(HH);
[I,J]=meshgrid([1:n],[1:m]);

IN=find(isnan(U1));
U1(IN)=0;
IN=find(isnan(V1));
V1(IN)=0;
IN=find(isnan(U2));
U2(IN)=0;
IN=find(isnan(V2));
V2(IN)=0;

% --------------------------------------------
% Runge-Kutta method
% Compute distance travelled by a particle
% and convert it into index space
% It is easier to work in index space because 
% of the 2D interpolation, index space is regular
% while Cartesian (km) are irregular
% much more difficult for interpolation
% --------------------------------------------
AX=[];
AY=[];
NP=length(LX);   % # of particles
for ip=1:NP
  ii=LX(ip);
  jj=LY(ip);
%keyboard
  % if outside domain - stop 
  if ii<=1 | ii>=n | jj<=1 | jj>=m
    AX(ip)=-100;
    AY(ip)=-100;
    continue;
  end

% --------------------------------------------
% k1x, k1y - velocity components at time level n, orig. location.
% --------------------------------------------
  u1=interp2(I,J,U1,ii,jj);
  v1=interp2(I,J,V1,ii,jj);
  k1x=u1;
  k1y=v1;
%keyboard
% --------------------------------------------
% Find k2x, k2y
% Interpolate U,V into point (t+0.5*dt,x+0.5*dt*k1x,y+0.5*dt*k1y)
% --------------------------------------------
  x1=k1x*0.5*dt;
  y1=k1y*0.5*dt;
% Translate distances x1 and y1 into index space:
  dxa=interp2(I,J,DX,ii,jj);
  ni=x1/dxa;
  dya=interp2(I,J,DY,ii,jj);
  nj=y1/dya;

  if1=ii+ni;
  jf1=jj+nj;

% if outside domain - stop 
  if if1<1 | if1>n | jf1<1 | jf1>m
    AX(ip)=-100;
    AY(ip)=-100;
    continue;
  end

% Time level n and n+1
% Interpolate in space:
  k2x_n1=interp2(I,J,U1,if1,jf1);
  k2y_n1=interp2(I,J,V1,if1,jf1);
  k2x_n2=interp2(I,J,U2,if1,jf1);  % time level n+1
  k2y_n2=interp2(I,J,V2,if1,jf1);

% Interpolate in time:
  k2x=0.5*(k2x_n1+k2x_n2);
  k2y=0.5*(k2y_n1+k2y_n2);

% --------------------------------------------
% Find k3x, k3y
% Time level 0.5*dt
% Interpolate U,V into point (t+0.5*dt,x+0.5*dt*k2x,y+0.5*dt*k2y)
% --------------------------------------------
  x2=k2x*0.5*dt;
  y2=k2y*0.5*dt;

% Translate distances x2 and y2 into index space:
  dxa=interp2(I,J,DX,if1,jf1);
  ni=x2/dxa;
  dya=interp2(I,J,DY,if1,jf1);
  nj=y2/dya;

  if2=ii+ni;
  jf2=jj+nj;

% if outside domain - stop 
  if if2<1 | if2>n | jf2<1 | jf2>m
    AX(ip)=-100;
    AY(ip)=-100;
    continue;
  end

% Interpolate in space:
  k3x_n1=interp2(I,J,U1,if2,jf2);
  k3y_n1=interp2(I,J,V1,if2,jf2);
  k3x_n2=interp2(I,J,U2,if2,jf2);  % time level n+1
  k3y_n2=interp2(I,J,V2,if2,jf2);

% Interpolate in time:
  k3x=0.5*(k3x_n1+k3x_n2);
  k3y=0.5*(k3y_n1+k3y_n2);

%-----------------------
% Find k4x, k4y
% Time level n+1
% ----------------------
  x3=k3x*dt;
  y3=k3y*dt;

  dxa=interp2(I,J,DX,if2,jf2);
  ni=x3/dxa;
  dya=interp2(I,J,DY,if2,jf2);
  nj=y3/dya;

  if3=ii+ni;
  jf3=jj+nj;

% if outside domain - stop 
  if if3<1 | if3>n | jf3<1 | jf3>m
    AX(ip)=-100;
    AY(ip)=-100;
    continue;
  end

% Interpolate in space:
  k4x=interp2(I,J,U2,if3,jf3);  % time level n+1
  k4y=interp2(I,J,V2,if3,jf3);

% Distance travelled by a particle:
  xnew=dt/6*(k1x+2*k2x+2*k3x+k4x);
  ynew=dt/6*(k1y+2*k2y+2*k3y+k4y);

% Get indices of the new location:
  dxa=interp2(I,J,DX,if3,jf3);
  ni=xnew/dxa;
  dya=interp2(I,J,DY,if3,jf3);
  nj=ynew/dya;

  inew=ii+ni;
  jnew=jj+nj;
  if inew<1 | inew>n | jnew<1 | jnew>m
    AX(ip)=-100;
    AY(ip)=-100;
    continue;
  end

% Add diffusion:
% Random walk 
% Gaussian
% xt=x(t-1)+u*dt+sqrt(2D*dt)*Xr
% yt=y(t-1)+v*dt+sqrt(2D*dt)*Yr
% Xr,Yr - random numbr Gauss. (0,1)
% D - diffusion, dt - time scale
% usually Dx>> Dy, if v=0
%
% In HYCOM, turbulent viscosity
% is given by (see HYCOM manual):
%  nu = max{ud*dltX, lmbd(div^2+vort^2)^1/2*dX^2}
%  ud is "diffusion velocity" that is nu/dltX
%  for momentum, ud is 2cm/s
% for 1/12dgr HYCOM, nu ~ 80 m2/s
%    nu = 80;  % m2/s
  sgm = sqrt(2*nu*dt);
  xdiff = sgm*randn(1);
  ydiff = sgm*randn(1);
  ni_diff=xdiff/dxa;
  nj_diff=ydiff/dya;
  
  iapr=inew+ni_diff;
  japr=jnew+nj_diff;
  if iapr<1 | iapr>n | japr<1 | japr>m
    iapr=inew;
    japr=jnew;
  end
% Check if particle is on shore
% if it is - move it away from coast:
%keyboard
  if isnan(japr) | isnan(iapr) |...
	iapr<1 | japr<1
    inew = -100;
    jnew = -100;
  else
    
  h0=HH(round(japr),round(iapr));
  if h0>-0.1;
%    keyboard
    di=sqrt((iapr-ii).^2+(japr-jj).^2); % prtcl displacement
    D=abs(sqrt((I-ii).^2+(J-jj).^2)-di);
    D(HH>=0)=nan;
%    D(D<di)=nan;
    [jdp,idp]=find(D==nanmin(nanmin(D)));
    inew=idp(1);
    jnew=jdp(1);
  else
    inew=iapr;
    jnew=japr;
  end
  end
  
%  keyboard
  
  AX(ip)=inew;
  AY(ip)=jnew;
end;  % for ip
%

return
