function [Xc,Yc,flg] = sub_combineLCLCE_EddyN(X,Y,LCE,irc,Neddy,varargin);
% Modified algorithm for combining LC and LCE number Neddy
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

flg=0;   % flag =1 if LCE has been added
Xc=X;
Yc=Y;
% No Eddies, keep LC contour only
if Neddy == 0, 
  return
end

% LCE is str array with LCEs
% limit to 1 LCEs to be combined - most eastward, LCE1>0 
nlce = length(LCE);

if ~isempty(lonW)
  if lonW>180; lonW=lonW-360; end;
end

if Neddy>nlce
  fprintf('sub_combineLCLCE: WARN: Eddy # %i > # of LCEs =%i \n',Neddy,nlce);
  return
end

%keyboard

ie = Neddy;
nee = length(LCE(ie).XY);
if irc>nee,
%	if ie==1  % # 
%	  fprintf('sub_combineLCLCE: requested time record outside record length %i %i\n',...
%			irc, nee);
%	end
  return
end

xe=LCE(ie).XY(irc).X;
ye=LCE(ie).XY(irc).Y;
if isempty(xe), % No LCE for this time record
  return;
end

% Find easternmost point of the LCE contour
eastX=max(xe);

if ~isempty(lonW) & eastX<lonW  % contour is west of cutoff
  return
end
  
Xc=[X;xe];
Yc=[Y;ye];
flg = 1;

return
