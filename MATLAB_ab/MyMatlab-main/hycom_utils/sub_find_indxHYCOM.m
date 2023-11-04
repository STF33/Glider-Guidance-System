function [IC,JC] = sub_find_indxHYCOM(XT,YT,xx,yy);
% function [IC,JC] = sub_find_indxHYCOM(XT,YT,xx,yy);
% For array of coordinates (xx,yy) Find corresponding indices 
% on HYCOM grid with coord. XT,YT
% use nearest point
np=length(xx);
for ii=1:np
  if mod(ii,100)==0
    fprintf('X,Y -> HYCOM indices: %6.2f%% done ...\n',ii/np*100);
  end
  x0=xx(ii);
  y0=yy(ii);
%  D=distance_spheric_coord(XT,YT,x0,y0);
  D=distance_spheric_coord(YT,XT,y0,x0);
  [j0,i0]=find(D==min(min(D)));
  IC(ii,1)=i0;
  JC(ii,1)=j0;
end


return
