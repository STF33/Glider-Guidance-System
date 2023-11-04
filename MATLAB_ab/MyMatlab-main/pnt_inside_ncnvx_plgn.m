         function [insd,nintrs]=pnt_inside_ncnvx_plgn(xx,yy,xp0,yp0,ray_drct)
% ==================================================
% Check if the point is inside:
% The easiest way is convex polygon and check if the point lies on one side 
% when going around the contour (see orientation.m)
% This does not work if the polygon is not convex
% Then Check # of crossings of a ray from the point with the contour lines
% If it is even (or zero) - point is outside the polygon, otherwise - inside
% The ray emitting from the point is parallel to X axis
% Special care should be taken to make sure the point is not one of
% the vertices and does not lay on the edge
% Input: xx,yy - coordinates of polygon vertices
%        xp0, yp0 - coordinates of the point
%
% Output: inds = 0 - the point is outside the polygon
%              = 1 - the point is inside the polygon
%         nintrs - # of crossings of the ray with the polygon edges
% ======================================================
  yintrs=yp0;   % ray is parallel to X
  nintrs=0;
  for mi=1:length(xx)-1
    aa1=xx(mi);
    aa2=xx(mi+1);
    bb1=yy(mi);
    bb2=yy(mi+1);
    if (aa2==aa1), aa2=aa1+1e-6; end;
    alf=(bb2-bb1)/(aa2-aa1);
    if (alf==0), alf=1e-6; end;
% Make sure that the point is not on the edge:
    ornt=orientation(aa1,bb1,aa2,bb2,xp0,yp0);
    if (ornt==0),
      disp('pnt_inside_ncnvx_plgn.m: '); 
      error('Point must not be on the edge of the polygon');
    end;
    xintrs=1/alf*(yintrs-bb1)+aa1;
    if (xintrs>=min(aa1,aa2) & xintrs<max(aa1,aa2) & xintrs>xp0), 
      nintrs=nintrs+1; 
    end;
 end;   % for mi
 if (floor(nintrs/2)*2==nintrs | nintrs==0),
   insd=0;
 else
   insd=1;
 end;









