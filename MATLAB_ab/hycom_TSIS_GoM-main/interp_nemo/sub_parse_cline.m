function aa2 = sub_parse_cline(aa,zz0,minv,maxv);
% Create c-line string from *b file
% for header in climat file Z-level HYCOM

strD = aa;
aa2 = sprintf('%s%8.1f%16.7e%16.7e',aa,abs(zz0),minv,maxv);

return

