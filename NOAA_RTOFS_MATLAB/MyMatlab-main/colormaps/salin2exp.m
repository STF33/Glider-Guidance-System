function AA = salin2exp(s0,S,flg);
% For plotting salinity - expand the 34 -35 range
% by transforming S into exp((s0-S).^(-1))
% flg = 1 - transform S -> exp()
% flg = -1 back conversion exp -> S
% Dmitry Dukhovskoy, COAPS FSU
%
if s0<=max(max(S)),
  fprintf('  !!!!!!!!!!!!!!!!!!!!! \n');
  warning('salin2exp: s0 should be > max S, s0=%6.4f, maxS=%6.4f',s0,max(max(S)));
  fprintf('  !!!!!!!!!!!!!!!!!!!!! \n');
end

if flg>0
  AA=exp((s0-S).^(-1));
else
  AA=s0-(log(S)).^-1;
end


return
