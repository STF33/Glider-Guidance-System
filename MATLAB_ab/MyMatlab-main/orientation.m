         function D = orientation(ax,ay,bx,by,cx,cy);
% ----------------------------------
%         function D = orientation(ax,ay,bx,by,cx,cy);
% Given vector ab by two points: A(ax,ay) - start point of vector
%  and B(bx,by) - end point, such that positive vector in X means
% bx-ax > 0
% and a point C(cx,cy)
% Does C lie on, to the left, to the right of ab?
% To find this out, find determinant of 
% | ax  ay  1 |
% | bx  by  1 | = | ax-cx  ay-cy |
% | cx  cy  1 |   | bx-cx  by-cy |
%
% From website of J. Shewchuk (author of triangle meshgrid)
%
% Output: determinate D
% If point ay < by, then D<0 means point C is to the right
%                        D>0 - point C is to the left
%                        D=0 - point C is on the line AB
% -------------------------------------------

  D = (ax-cx)*(by-cy)-(ay-cy)*(bx-cx);
