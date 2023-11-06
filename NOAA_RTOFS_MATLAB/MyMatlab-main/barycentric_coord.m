        function [lmb1,lmb2,lmb3] = barycentric_coord(XYT,x0,y0);

% Find barycentric coordinates of a point for a given triangle with 
% vertices: xm1,ym1; xm2,ym2; xm3,ym3
%! lambdas should be >0 and <1
% Note that if any of lambdas (lmb1,...) is > 1, the point is outside 
% the triangle
% if one of the lmbds ==0 or ==1 - point lies on the edge or corner of a triangle
%
% given: XYT=(xm1,xm2; ...) vertices, coordinates
%        x0,y0  - coordinates of a point
% 
% Output: lmb1, lmb2, lmb3 - barycentric coordinates 
%

xm1=XYT(1,1);
ym1=XYT(1,2);
xm2=XYT(2,1);
ym2=XYT(2,2);
xm3=XYT(3,1);
ym3=XYT(3,2);

r=[x0;y0];
r3=[xm3;ym3];
T=[xm1-xm3, xm2-xm3; ym1-ym3, ym2-ym3];
% 
% Check if matrix is singular, meaning
% that 3 pnts lie on 1 line
if abs(det(T))<1e-9
  lmb1=-100;
  lmb2=-100;
  lmb3=-100;
else
  lmb12=inv(T)*(r-r3);
  lmb1=lmb12(1);
  lmb2=lmb12(2);
  lmb3=1-lmb1-lmb2;
end

