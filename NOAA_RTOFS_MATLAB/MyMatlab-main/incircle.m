  function   k = incircle(ax,ay,bx,by,cx,cy,dx,dy);
% Algoritm of incircle test:
% does point D(dx,dy) lie inside (k>0), or
% outside (k<0) of ABC?
% Based on Shewchuk, 
% http://www.cs.cmu.edu/~quake/robust.html

DD=[ax,ay,(ax^2+ay^2),1;...
    bx,by,(bx^2+by^2),1;...
    cx,cy,(cx^2+cy^2),1;...
    dx,dy,(dx^2+dy^2),1];

  k = det(DD);
