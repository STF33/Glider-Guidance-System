  function sub_check_B(nlev,SD,td0,minv,maxv);
% Check read in *b file with
% output information
% Layer #, targ. dens. (td0) should match
%  SD - values from *b string, at level nlev
%       parsed from a string read from *b file
%  that typically has format:
% thknss   =          1      15.00  2 20.250   9.8060000E+03   7.5407562E+04
% SD - values after "="
%
dmn = abs(minv-SD(end-1));  
dmx = abs(maxv-SD(end));
if dmn>0.01| dmx>0.01
  error('*** ERR:  min/max values do not agree...');
end

dd=abs(SD(end-2)-td0);
if dd>1e-3
  fprintf('Target Dens for k=%i is %f\n',nlev,td0);
  fprintf('From *b file it is %f\n',SD(end-2));
  error('*** STOPPING');
end


return