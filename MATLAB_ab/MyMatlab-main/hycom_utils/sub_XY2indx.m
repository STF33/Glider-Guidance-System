function IJ = sub_XY2indx(XY,LON,LAT);
% IJ = sub_XY2indx(XY,LON,LAT);
% XY - pair of X and Y coordinates
% Find corresponding indices
np = size(XY,1);
rmd = round(np/10);
for ip=1:np
  if rmd>100 & mod(ip,rmd)==0
    fprintf('sub_XY2indx: %5.1f%% done ...\n',ip/np*100);
  end
  x=XY(ip,1);
  y=XY(ip,2);
  
  dd = distance_spheric_coord(y,x,LAT,LON);
  [j,i] = find(dd==min(min(dd)),1);
  IJ(ip,1)=i;
  IJ(ip,2)=j;
end



return
