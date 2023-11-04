function [DX,DY]=sub_dx_dy(LON,LAT);
% [DX,DY]=sub_dx_dy(LON,LAT);
% for given LON & LAT arrays
% calculated dx and dy for each grid cell
%
%
[mm,nn]=size(LON);
DX=zeros(mm,nn);
DY=zeros(mm,nn);
for i=1:nn-1
  dx=distance_spheric_coord(LAT(:,i),LON(:,i),LAT(:,i+1),LON(:,i+1));
  DX(:,i)=dx;
end
DX(:,nn)=dx;
for j=1:mm-1
  dy=distance_spheric_coord(LAT(j,:),LON(j,:),LAT(j+1,:),LON(j+1,:));
  DY(j,:)=dy;
end
DY(mm,:)=dy;

return
