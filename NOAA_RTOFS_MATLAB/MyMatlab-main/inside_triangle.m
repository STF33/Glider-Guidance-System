     function k = inside_triangle(ax,ay,bx,by,cx,cy,x0,y0);
% function k = inside_triangle(ax,ay,bx,by,cx,cy,x0,y0);
% Deterimnes if point (x0,y0) lies inside (k > 0), outside (k < 0)
% or on one of the sides (k = 0) of triangle with nodes A(ax,ay), B(bx,by),
% and C(cx,cy)
% I use test: for two lines starting at one of the nodes the point
% should lie on different sides if it is inside the triangle
% The test should be run for any two different nodes
% To determine the orientation of the node relative to a line
% I use algorythm of Shewchuk
  A1 = orientation(ax,ay,bx,by,x0,y0);     
  A2 = orientation(ax,ay,cx,cy,x0,y0);
  if (abs(A1)<1e-8) | (abs(A2)<1e-8)
    sA=0;
  else
    sA = A1/abs(A1)*A2/abs(A2);
  end;

  B1 = orientation(bx,by,ax,ay,x0,y0);     
  B2 = orientation(bx,by,cx,cy,x0,y0);
  if (abs(B1)<1e-8) | (abs(B2)<1e-8)
    sB=0;
  else
    sB = B1/abs(B1)*B2/abs(B2);
  end;
  
  if (abs(sA)<1e-8) | (abs(sB)<1e-8)
    k=0;                   % on a side
  elseif (sA<0) & (sB<0)
    k=1;                   % inside 
  else
    k=-1;                  % outside
  end;
  
