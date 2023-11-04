% Find indices for subregions in the GOM
% 501 GOMu reanalysis
function GOMR = sub_find_501GOMindx(GOM,LON,LAT,lonR,latR);

nvrt = size(GOM,1);
for ii=1:nvrt
  i0 = GOM(ii,1);
  j0 = GOM(ii,2);
  x0 = LON(j0,i0);
  y0 = LAT(j0,i0);

% Assumed that the domain does not cross 0 meridian
  ir1 = max(find(lonR<=x0));
  jr1 = max(find(latR<=y0));

  if y0 < min(latR),
    jr1 = 1;
  end
  if x0 < min(lonR),
    ir1 = 1;
  end
  if y0 > max(latR),
    jr1 = length(latR);
  end
 
  GOMR(ii,1) = ir1;
  GOMR(ii,2) = jr1;

end

return 


