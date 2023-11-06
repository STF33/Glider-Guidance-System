   function rmsd = sub_rmsd(P,Q);
%
% Root mean square deviation (error) is
% clacluated for 2 contours Q and P
%
%
% Need to generalize code when 
% P and Q are of different sizes,
% then exgtra points are added along
% segments of the contour with less elements  
%
% Points have to be ordered properly, i.e.
% P(1,1) is assumed to correspond to Q(1,1)
% etc.
%
%
%

[a1,a2]=size(P);
if a2>a1, 
  P=P';
end;
[a1,a2]=size(Q);
if a2>a1, 
  Q=Q';
end;

np=size(P,1);
nq=size(Q,1);

d=np-nq;
addP=0;
addQ=0;
if d~=0

  error('Need to finish this code for Q~=P size');
  
  addC=abs(d);
  if d>0 % add to Q
    add1sgm=addC/nq; %how many pnts to add per 1 segm
    dmm=Q;
    D0=P;
  else
    add1sgm=addC/np; %how many pnts to add per 1 segm
  end

  nadded=0;
  dadd=0;
  for j=1:length(D0)-1;
    if dadd>1
      npnts=floor(dadd);
      x1=Q(j,1);
      y1=Q(j,2);
      x2=Q(j+1,1);
      y2=Q(j+1,2);
      ll=sqrt((x1-x2).^2+(y1-y2).^2);
      lss=ll/(npnts+1);
%
% Not finished here ...
      
      dadd=dadd-npnts;
    end
    
    dadd=dadd+add1sgm;
  end
  

end

D=(P(:,1)-Q(:,1)).^2+(P(:,2)-Q(:,2)).^2;
rmsd=sqrt(sum(D)./np);

