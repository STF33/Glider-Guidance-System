function [Xc,Yc] = sub_combineLCLCE_v2(X,Y,LCE,irc,LCE1);
% Modified algorithm for combining LC and LCE
% The closes LCE is chosen 
% the closest point on LCE and LC is found 
%
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
    npp  = length(xe);
    minD  = 1e20;
    iLCE  = [];  % LCE
    iLC   = [];  % index on LC
    for ipp=1:npp
      d=sqrt((xe(ipp)-X).^2+(ye(ipp)-Y).^2);
      md=min(d);
      minD=min([md,minD]);
      if minD==md
        iLCE=ipp;
        iLC=find(d==md);
      end
    end
    k0=find(D==min(D),1);
    clsX(ie,1)=xe(iLCE);
    clsY(ie,1)=ye(iLCE);
    clsX(ie,2)=X(iLC);
    clsY(ie,2)=Y(iLC); 
    dstLCE(ie)=minD;  % closest distance
  end
%keyboard    

  if LCE1~=1 
    X=[X;xe];
    Y=[Y;ye];
  end
end

%keyboard

if isempty(eastX),  % no LCE 
  Xc=X;
  Yc=Y;
  return; 
end; 


% If >1 LCEs, keep the closest 
% to the LC, i.e the most eastward
if LCE1==1 & ~isempty(dstLCE)
  iE = find(dstLCE==min(dstLCE));
		xe=LCE(iE).XY(irc).X;
		ye=LCE(iE).XY(irc).Y;
		X=[X;xe];
		Y=[Y;ye];
end  
  
Xc=X;
Yc=Y;



return
