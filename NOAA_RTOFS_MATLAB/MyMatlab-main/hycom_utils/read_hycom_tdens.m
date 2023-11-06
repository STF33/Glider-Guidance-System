    function TDENS = read_hycom_tdens(finb);
% 
% reads target densities from *.b file
% l - # of target densities
%
f1 = fopen(finb,'r');  % read I,J from *.b
for nl=1:1e6
  aa=fgetl(f1);
%  disp(aa);
  if ~ischar(aa); error('Target Dens: Couldnot locate fields'); end
  I= strfind(aa,'u-vel.');
  if ~isempty(I); break; end;
end

disp(['Reading target densities ', finb]);
%
a=1;
kk=0;
cc=0;
while a
  aa=fgetl(f1);
  if ~ischar(aa); break; end

  [ss1, ss2, tstp, mday, k, tdns, amm, axx] = strread(aa,'%s%s%d%f%d%f%f%f');
  if k>kk
    kk=k;
    cc=cc+1;
    TDENS(cc,1)=tdns;
    K(cc,1)=k;
  end
%  if (k>37); keyboard; end;

end
fclose(f1);
fprintf('Read in target densities = %i \n',cc);
