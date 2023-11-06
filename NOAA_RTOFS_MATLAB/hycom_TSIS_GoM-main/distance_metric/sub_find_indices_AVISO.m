% Find indices of the GoM in AVISO
function [ia1,ia2,ja1,ja2] = sub_find_indices_AVISO(flaviso,xt1,xt2,yt1,yt2);

fprintf('Reading AVISO grid %s\n',flaviso);
lonA = nc_varget(flaviso,'longitude');
latA = nc_varget(flaviso,'latitude');

I = find(lonA>180.);
lonA(I) = lonA(I)-360.;

% Assumed that the domain does not cross 0 meridian
ia1 = max(find(lonA<=xt1));
ia2 = max(find(lonA<=xt2))+1;
ja1 = max(find(latA<=yt1));
ja2 = max(find(latA<=yt2))+1;

fprintf('AVISO /HYCOM Min/max lon = %5.2f/%5.2f  %5.2f/%5.2f\n',...
        lonA(ia1),lonA(ia2),xt1,xt2);
fprintf('AVISO /HYCOM Min/max lat = %5.2f/%5.2f  %5.2f/%5.2f\n',...
        latA(ja1),latA(ja2),yt1,yt2);

if ia2<ia1
  error(' second lon index (%i) < 1st (%i), 0-meridian problem?',ia2,ia1);
end

return
