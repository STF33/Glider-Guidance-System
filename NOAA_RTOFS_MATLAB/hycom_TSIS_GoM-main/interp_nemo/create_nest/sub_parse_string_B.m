function aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,varargin);
% Update string from *b file
% Replace min max values

nV = length(varargin);
hycom_date=0;
for k=1:nV
  scl=varargin{k};
  scl=lower(scl);
  if strncmp(scl,'hycom_date',8)
    hycom_date=varargin{k+1};  % hycom model date
  end
end
%keyboard
strD = aa;
I=strfind(aa,'=');
bb=aa(I+1:end);
S = sscanf(bb,'%f');
if hycom_date>0
  S(2)=hycom_date;
end

if isempty(minv) | isempty(maxv)
  minv = S(5);
  maxv = S(6);
end

% String: info, t.step, mean day, vert. layer#, t. dens., min va., max val
if ~isempty(nlev) & ~isempty(td0)
  aa2 = sprintf('%s%11i%11.2f%3i%7.3f%16.7e%16.7e',...
	       aa(1:I),S(1),S(2),nlev,td0,minv,maxv);
else
  aa2 = sprintf('%s%11i%11.2f%3i%7.3f%16.7e%16.7e',...
	       aa(1:I),S(1),S(2),S(3),S(4),minv,maxv);
end

return
