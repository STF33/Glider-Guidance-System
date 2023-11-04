% Define vertical section to plot T/S profiles
% Eastern GoM 
function SCT = sub_GoM_xsct(LON,LAT,HH);
SCT.XY=[-86.44657, 20.13259;...
       -84.89643, 23.36478; ...
       -90.0745, 27.553];

nsgm = size(SCT.XY,1);
for isgm=1:nsgm
  XY = SCT.XY(isgm,:);
  IJ = sub_XY2indx(XY,LON,LAT);
  SCT.IJ(isgm,:)=IJ;
end

IJs=SCT.IJ;
nij=size(IJs,1);
IIs=[];
JJs=[];
for ii=1:nij-1
  i1=IJs(ii,1);
  i2=IJs(ii+1,1);
  j1=IJs(ii,2);
  j2=IJs(ii+1,2);
  [I,J]=sub_xsct_indx(i1,j1,i2,j2);
  if size(I,1)==1;
    I=I';
    J=J';
  end

  if isempty(IIs)
    IIs=[IIs;I];
    JJs=[JJs;J];
  else
    IIs=[IIs;I(2:end)];
    JJs=[JJs;J(2:end)];
  end
end;

SCT.I=IIs;
SCT.J=JJs;


IIs=SCT.I;
JJs=SCT.J;
nsg=length(IIs);
clear XX YY
for ii=1:nsg
  i0=IIs(ii);
  j0=JJs(ii);
  XX(ii)=LON(j0,i0);
  YY(ii)=LAT(j0,i0);
  Hbtm(ii)=HH(j0,i0);
end

INDs = sub2ind(size(HH),JJs,IIs);

SCT.Indx = INDs;
SCT.long = XX';
SCT.latd = YY';
SCT.Hbtm = Hbtm';

% Compute segments' lengths:
lii = length(XX);
DST=zeros(lii,1);
for ii=1:lii-1
  x1=SCT.long(ii);
  y1=SCT.latd(ii);
  x2=SCT.long(ii+1);
  y2=SCT.latd(ii+1);

  dd=distance_spheric_coord(y1,x1,y2,x2)*1e-3; % km
  DST(ii+1)=dd;
end
SCT.dX = DST;

return

