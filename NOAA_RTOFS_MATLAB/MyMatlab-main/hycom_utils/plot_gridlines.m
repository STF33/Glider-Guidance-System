  function plot_gridlines(dlmb,dphi, nf, clr, elon, alat)
% 
% Draw long/lat on the map
% will need elon and alat
%
% Input; 
%   dlmb - delta lambda, longitudes lines 
%   dphi - delta latitudes circles
%   nf   - figure # where to plot
%   clr - color of the coordinate lines
I=find(elon>180);
elon(I)=elon(I)-360;
elon(alat>=89.)=nan;
el1=elon;
el1(elon>178.)=nan;
el1(elon<-178.)=nan;
contour(el1,[-180:dlmb:180],'Color',clr,'linewidth',1);
el1=elon;
I=find(el1<0);
el1(I)=el1(I)+360;
el1(el1<50)=nan;
el1(el1>300)=nan;
contour(el1,[180 180],'Color',clr,'linewidth',1);

[rr1,gg1]=contour(alat,[40:dphi:90],'Color',clr);

return



