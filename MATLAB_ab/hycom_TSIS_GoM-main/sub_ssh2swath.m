% get SSH values in satellite swaths coordiantes
% use nearest coord from SSH grid
% within the domain: lon1/lon2 - lat1/lat2
%
function AA = sub_ssh2swath(X,Y,SSH,LON,LAT,lon1,lon2,lat1,lat2);

nx = length(X);
icc=0;
for ix=1:nx
  x0=X(ix);
  y0=Y(ix);

  if x0<lon1 | x0>lon2 | y0<lat1 | y0>lat2
    continue;
  end

  D=distance_spheric_coord(LAT,LON,y0,x0);
  [j0,i0] = find(D==min(min(D)));

% Check if same grid point is already selected:
  if icc>0
    II = find(AA.Iindx==i0 & AA.Jindx==j0);
    if ~isempty(II); continue; end;
  end

  icc=icc+1;
  AA.X(icc)      = x0;
  AA.Y(icc)      = y0;
  AA.Iindx(icc)  = i0;
  AA.Jindx(icc)  = j0;
  AA.lonH(icc)   = LON(j0,i0);
  AA.latH(icc)   = LAT(j0,i0);
  AA.ssh(icc)    = SSH(j0,i0);
end

