function [Xc,Yc,flg] = sub_combineLCLCE_EComb(X,Y,LCE,irc,EComb,lcomb,varargin);
% Combine LC and LCEs using EComb struct array (sub_lce_combm)
% with all possible combinations of LCEs - pick selected combination lcomb
% 
lonW=[];
nV = length(varargin);
if nV > 0
  for k=1:nV
    vfld = varargin{k};
    if strmatch(vfld,'lonW')  % longitude westward Eddies are discarded 
      lonW = varargin{k+1};
    end
  end
end
if ~isempty(lonW)
  if lonW>180; lonW=lonW-360; end;
end


flg=0;   % flag > 1 if LCE has been added
Xc=X;
Yc=Y;

icc=0;
lce = [];
ll=0;
while icc < lcomb
  ll=ll+1;
  aa=EComb(ll).LCE;
  id1 = size(aa,1);
  id2 = size(aa,2);
  for ilce=1:id1
    icc=icc+1;
    if icc == lcomb; 
      lce = aa(ilce,:);
      break; 
    end;
  end
end

% No Eddies, keep LC contour only
if isempty(lce), 
  fprintf('!! sub_combine: Comb # mismatch ihc=%i, lcomb=%i\n\n',irc,lcomb);
  error('*ERR: Check lcomb, EComb ')
end

if lce == 0; return; end;  % no LCEs

% LCE is str array with LCEs
for ikk=1:length(lce)
  ie = lce(ikk);
  nee = length(LCE(ie).XY);
  if irc>nee,
    flg = -1;
    return
  end

  xe=LCE(ie).XY(irc).X;
  ye=LCE(ie).XY(irc).Y;
  if isempty(xe), % No LCE for this time record
    flg = -1;
    return;
  end

% Find easternmost point of the LCE contour
  eastX=max(xe);

  if ~isempty(lonW) & eastX<lonW  % contour is west of cutoff
    continue;
  end

  flg = flg+1;
  Xc=[Xc;xe];
  Yc=[Yc;ye];
end

%flg = ikk;
if flg < length(lce)
  flg = -ikk;  % some LCEs were not added (west of lonW) 
end

return
