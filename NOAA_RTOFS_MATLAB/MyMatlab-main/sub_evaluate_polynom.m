function   [gs,xxs] = sub_evaluate_polynom(Xis,Yis,Fdf,xx,dn);
% Evaluate polynomial at points Xis
% using precomputed divided differences Fdf (Smoktunowicz)
% Xis, Yis - x,y data points in the segment
% xx - points where polynom is evaluated
% dn - polynomial degree in the segment
%
gs = [];
xxs = [];
x0s = Xis(1);
xNs = Xis(end);
ix1 = min(find(xx>=x0s & xx<xNs));
ix2 = max(find(xx<=xNs & xx>x0s));
if isempty(ix1) | isempty(ix2)
  fprintf('Segment %i contains no data points\n',isgm);
  return
end
xxs = xx(ix1:ix2); % subsample xx points for this segment
xxs = xxs(:);

%fprintf('ix1=%i ix2=%i\n',ix1,ix2);

dmm = Fdf(1); % y0
dff = 1;

for ii=2:dn+1
  dd = xxs-Xis(ii-1); % (x-x0)*(x-x1)...
  dff = dff.*dd;
  dmm = dmm+Fdf(ii).*dff;
end
gs = dmm;
gs=gs(:);

return
  