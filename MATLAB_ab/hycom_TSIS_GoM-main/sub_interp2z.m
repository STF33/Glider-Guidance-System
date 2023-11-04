% Interpolate from HYCOM hybrid grid 
% onto vertical z-grid with fixed depth layers
function Vi = sub_interp2z(V,ZZ,ZZf);

nint = length(ZZf);
nx = size(V,2);
nz = size(ZZ,1);
Vi = zeros(nint,nx);
% Add missing V for bottom interface
V(end+1,:) = 0;
for ix=1:nx
% fill in bottom:
  zz = ZZ(:,ix);
  vv = V(:,ix);
  Ib = min(find(isnan(zz)));
  if nanmin(zz)==0 | isnan(vv(1)); continue; end;
  for iz=Ib:nz
    zz(iz)=zz(iz-1)+1.e-8;
    vv(iz)=0.;
  end

% for interpolation - add super deep point
  zz(end+1)=-8000;
  vv(end+1)=0.;
  vvi = interp1(zz,vv,ZZf);

  Vi(:,ix)=vvi;
end
  
 

return
