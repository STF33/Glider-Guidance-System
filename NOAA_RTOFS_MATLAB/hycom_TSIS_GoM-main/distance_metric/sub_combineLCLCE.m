function [Xc,Yc] = sub_combineLCLCE(X,Y,LCE,lnW,irc,LCE1);
% Combine LC and LCE contours east of cutoff lon lnW
% LCE is str array with LCEs
% limit to 1 LCEs to be combined - most eastward, LCE1>0 
nlce = length(LCE);

maxX = [];
mnX = [];

if nlce==0
  Xc=X;
  Yc=Y;
  return
end

if lnW>180, lnW=lnW-360; end
if max(X)>0; error('sub_combineLCLCE: Longitude in GoM should be <0'); end;


% Farthest point from YC
x0=-85;
y0=22;
Xmm=X;
Ymm=Y;
Xmm(Xmm>-84)=nan;  % chop off eastern part in Str. Florida
D=sqrt((Xmm-x0).^2+(Ymm-y0).^2);
j0=find(D==max(D));
xLC0=X(j0);
yLC0=Y(j0);
dstLCE=[];
clsX=[];
clsY=[];
eastX=[];
northY=[];
for ie=1:nlce
  nee = length(LCE(ie).XY);
		if irc>nee,
    if ie==1  % # of records for LC should be full lengths
    fprintf('sub_combineLCLCE: requested time record outside record length %i %i\n',...
      irc, nee);
    end
				continue;
		end

		xe=LCE(ie).XY(irc).X;
		ye=LCE(ie).XY(irc).Y;
		if isempty(xe),
				continue;
		end
  eastX(ie)=max(xe);
  northY(ie)=max(ye);
% Find closest point on the LCE eddy to the LC
%keyboard
  if LCE1==1 & ~isempty(xe)
    D=sqrt((xe-xLC0).^2+(ye-yLC0).^2);
    k0=find(D==min(D),1);
    clsX(ie)=xe(k0);
    clsY(ie)=ye(k0);
    dstLCE(ie)=D(k0);  % closest distance
  end
%keyboard    

  if LCE1~=1 & eastX>lnW
    X=[X;xe];
    Y=[Y;ye];
  end
end

%keyboard

% Chop the LC west of lnW
if isempty(eastX),  % no LCE 
%  II=find(X>=lnW);
%  Xc=X(II);
%  Yc=Y(II);
  Xc=X;
  Yc=Y;
  return; 
end; 


% If >1 LCEs, keep the closest 
% to the LC, i.e the most eastward
if LCE1==1 & ~isempty(dstLCE)
  iE = find(dstLCE==min(dstLCE));
%
% only if LCE is east of lnW
  if ~isempty(lnW) 
    if eastX(iE)>=lnW    
      xe=LCE(iE).XY(irc).X;
      ye=LCE(iE).XY(irc).Y;
      X=[X;xe];
      Y=[Y;ye];
    end
  end
end  
  

%II=find(X>=lnW);
%Xc=X(II);
%Yc=Y(II);
Xc=X;
Yc=Y;



return
