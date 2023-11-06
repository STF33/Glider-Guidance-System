    function [X,Y]=get_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
% ----------------------------------
% get_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
%
% draws arrow with filled "head"
% and returns coordinates of the arrow head
% This could be useful if other drawing 
% needs to be done with the vector
%
% Note: first X-coordinates x1, x2, then Y-coordinates, Y1,y2, ...
% Draws an arrow between pnt (x1,y1) and pnt (x2,y2)
% cf - scaling coefficient of the arrowhead, if [0 to 1)
%      arrow head is smaller than the vector
%      otherwise, the vector will be closed by the arrowhead
% beta - angle between vector and arrow head beams (degrees)
% v_col - color ([R G B]) % default color is black
% lwd - line width, default = 1
% 
 if (isempty(v_col)); v_col=[0,0,0]; end; % default color is black
 if (isempty(lwd)); lwd=1.; end; % default line width

  hold on
  uu=x2-x1;
  vv=y2-y1;
  sp=sqrt(uu.*uu+vv.*vv);	  
  alfa=atan2(uu,vv);		      % vector angle from Y
  beta=beta*pi/180;
  var=cf*sp;				          % scaling of the arrow head
  dX2=var.*sin(alfa-beta);		  % arrow head coordinates
  dX3=var.*sin(alfa+beta);		  % 
  dY2=var.*cos(alfa-beta);
  dY3=var.*cos(alfa+beta);
  dL=sqrt(dX2^2+dY2^2);
%  dL=sqrt(dX3^2+dY3^2);
%
% Length of the vector with the arrow-head
  Lv=sp+dL-0.02*dL; % to avoid gap btw stem and arrowhead
  un=uu/sp;
  vn=vv/sp;
  x0=x1+un*Lv; % scale to adjust for the arrowhead
  y0=y1+vn*Lv;
  
  ax2=x0-dX2;
  ax3=x0-dX3;
  ay2=y0-dY2;
  ay3=y0-dY3;
%keyboard
%  x2v=x1+(1-cf)*uu;
%  y2v=y1+(1-cf)*vv;
%  p1=plot([x1 x2v],[y1 y2v],'Color',v_col);     %vector
  p1=plot([x1 x2],[y1 y2],'Color',v_col);     %vector
%  p2=plot([x2 ax2],[y2 ay2],'Color',v_col);		% arrow head
%  p3=plot([x2 ax3],[y2 ay3],'Color',v_col);
 
  set(p1,'linewidth',lwd);
%  set([p2,p3],'linewidth',0.5);

    X=[x0,ax2,ax3,x0];
    Y=[y0,ay2,ay3,y0];
    H=fill(X,Y,v_col);
    set(H,'edgecolor',v_col);

%keyboard

%  hold off
return
