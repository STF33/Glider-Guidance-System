% Find LCE combinations given the number of LCEs
% For each number of eddies - given combination of LCEs + LC
%
function [EComb, icomb] = sub_lce_comb(Nlc2);

icomb = 0;
EComb = [];

% No eddies, only LC
ie = 1;
EComb(ie).Name = sprintf('LC + %i LCEs',ie-1);
EComb(ie).LCE = 0;
icomb = icomb+1;

if Nlc2 == 0;
  return
end

ie = ie+1;
EComb(ie).Name = sprintf('LC + %i LCEs',ie-1);
for ll=1:Nlc2
  EComb(ie).LCE = [1:Nlc2]'; % 1 eddy
  icomb = icomb+1;
end

if Nlc2 == 1; return; end

% Comb of 2 LCEs
ie = ie+1;
EComb(ie).Name = sprintf('LC + %i LCEs',ie-1);
icc=0;
for ii=1:Nlc2
  for jj=ii:Nlc2
    if ii ~= jj
      icc = icc+1;
      EComb(ie).LCE(icc,1:2) = [ii,jj];
      icomb = icomb+1;
    end
  end
end

if Nlc2 == 2; return; end;

% Comb of 3 LCEs
ie = ie+1;
EComb(ie).Name = sprintf('LC + %i LCEs',ie-1);
icc=0;
for ii=1:Nlc2
  for jj=ii:Nlc2
    for kk=jj:Nlc2
      if ii ~= jj & ii ~= jj & jj~= kk
        icc = icc+1;
        EComb(ie).LCE(icc,1:3) = [ii,jj,kk];
        icomb = icomb+1;
      end
    end
  end
end

if Nlc2 == 3; return; end;

% Comb of 4 LCEs
ie = ie+1;
EComb(ie).Name = sprintf('LC + %i LCEs',ie-1);
icc=0;
for ii=1:Nlc2
  for jj=ii:Nlc2
    for kk=jj:Nlc2
      for ll=kk:Nlc2
        if ii ~= jj & ii ~= kk & ii ~= ll & ...
           jj ~= kk & jj ~= ll & ...
           kk ~= ll
          icc = icc+1;
          EComb(ie).LCE(icc,1:4) = [ii,jj,kk,ll];
          icomb = icomb+1;
        end
      end
    end
  end
end

% No more than 4 LCEs are assumed

if Nlc2 > 4;
  fprintf('Nlc2=%i > Assumed max number of LCEs is 4\n',Nlc2);
  error('Need to modify sub_lce_comb.m');
end


return 

