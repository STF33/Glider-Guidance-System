function aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv);
% Update string from *b file
% Replace min max values

strD = aa;
I=strfind(aa,'=');
bb=aa(I+1:end);
S = sscanf(bb,'%f');
% String: info, t.step, mean day, vert. layer#, t. dens., min va., max val
if ~isempty(nlev) & ~isempty(td0)
  aa2 = sprintf('%s%11i%11.3f%3i%6.2f%16.7e%16.7e',...
	       aa(1:I),S(1),S(2),nlev,td0,minv,maxv);
else
%  aa2 = sprintf('%s%11i%11.2f%3i%7.3f%16.7e%16.7e',...
  aa2 = sprintf('%s%11i%11.3f%3i%6.2f%16.7e%16.7e',...
	       aa(1:I),S(1),S(2),S(3),S(4),minv,maxv);
end

return