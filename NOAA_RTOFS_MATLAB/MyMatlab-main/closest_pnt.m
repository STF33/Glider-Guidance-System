    function [i0,j0] = closest_pnt(X,Y,x0,y0);
%    function [i0,j0] = closest_pnt(X,Y,x0,y0);
%
%  function finds the closest point to (x0,y0) in the 1D or 2D array 
%  of X,Y coordinates 
%

dst=sqrt((X-x0).^2+(Y-y0).^2);
[j0,i0]=find(dst==min(min(dst)));



