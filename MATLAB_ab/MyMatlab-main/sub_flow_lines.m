function sub_flow_lines(U,V,X,Y,DX,DY,CF,HH);
% Draws flow lines based on 
% provided velocity field
% under steady-state assumption
% Input: U,V - velocities
%       X,Y - coordinates for plotting
%       DX,DY - grid distances in meters
%  CF - structured array with
%       figure properties and 
%       flow line specifications
%
%addpath /Users/dmitry/Documents/work/codes/MyMatlab;
%addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
%startup;


[mm,nn]=size(X);
% Steady-state: flow fields does not change:
U1=U;
U2=U;
V1=V;
V2=V;
%DX=X*0+500;
%DY=DX;
%HH=DX*0-100;
nu=0;
%[LX,LY]=meshgrid([1:4:nn],[1:4:mm]);
%[a1,a2]=size(LX);
%np=a1*a2;
%LX=reshape(LX,a1*a2,1);
%LY=reshape(LY,a1*a2,1);
% Locations where flow lines are plotted
% These are INDICES, 1D arrays
LX = CF.initial_flow_I;
LY = CF.initial_flow_J;
np=length(LX);
dt=800;
clear XP YP
XP(1,:)=LX;
YP(1,:)=LY;
dt = CF.Runge_Kutta_dT; %%%%%%%%% time stepping
iend = CF.Runge_Kutta_Nitir;

for ik=1:iend
  fprintf('ik=%i\n',ik);
  [AX,AY]=runge_kutta(U1,V1,dt,U2,V2, LX, LY, DX, DY,HH,nu);
  XP(ik+1,:)=AX;
  YP(ik+1,:)=AY;
  LX=AX;
  LY=AY;
%  fprintf('Max X=%10.2f, min X=%10.2f\n',max(XP),min(XP));
end

XP(XP<0)=nan;
YP(YP<0)=nan;

cf=CF.cf;
beta=CF.beta;
v_col=CF.v_col;
lwd=1.;

% Interpolate into coordinnates X,Y 
% from indeces
[II,JJ]=meshgrid([1:nn],[1:mm]);
XPi=interp2(X,XP,YP);
YPi=interp2(Y,XP,YP);

%keyboard
% Plot streamlines:
%%figure(2); clf;
%hold on;
for ik=1:np
  a=XPi(:,ik);
  b=YPi(:,ik);
  plot(a,b);
% Plot arrow-heads
  I=max(find(~isnan(a)));
  if I<2, continue; end;
  x1=a(I-1);
  x2=a(I);
  y1=b(I-1);
  y2=b(I);
% Scale to make the arrow  = 1
  dx=(x2-x1);
  dy=(y2-y1);
  L=sqrt(dx^2+dy^2);
  dx=dx/L;
  dy=dy/L;
  x1c=x1-dx;
  y1c=y1-dy;
  draw_arrow(x1c,x2,y1c,y2,cf,beta,v_col,lwd);
end;
axis('equal');
return


