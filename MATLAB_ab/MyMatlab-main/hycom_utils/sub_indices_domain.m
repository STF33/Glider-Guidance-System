function IJ = sub_indices_domain(XY,LON,LAT);
% function IJ = sub_indices_domain(XY,LON,LAT);
% COnvert Lon,Lat coordinates in XY(N,2) array
% into indices for given LON, LAT
n=size(XY,1);
fprintf('sub_indices_domain: Searching for indices, %i points ...\n',n);
for ik=1:n
  x0=XY(ik,1);
  y0=XY(ik,2);
  d=distance_spheric_coord(y0,x0,LAT,LON);
  [j,i]=find(d==min(min(d)),1);
  IJ(ik,1)=i;
  IJ(ik,2)=j;
end


return