  function [Uc,Vc] = sub_putUVonCgrid(UU,VV,HH);
% U and V in the nest files
% are on C-grid:
%     -------------
%     |           |
%     |           |
%     u    dP(i,j)  
%     |           |
%     |           |
%     ----- v ------
%
%  This is particularly important near the land and OBs
%
%  OBs: V are all 2^100 at South, East, North
%       U is 2^100 at East & North only
%
% Here, I simply shift U and V by half grid
% and make nans at the land boundaries
%
hg     = 2^100;

[ll,mm,nn]=size(UU);

Uc=UU*nan;
Vc=VV*nan;

for k=1:ll
  u=squeeze(UU(k,:,:));
  v=squeeze(VV(k,:,:));
  for i=2:nn-1
    h=HH(:,i-1);
    Ld=find(h>=0);
    u=squeeze(UU(k,:,i));
    if isnan(u(1))& h(1)<0, u(1)=u(2); end;
    u(end)=nan;
    u(Ld)=nan; 
    Uc(k,:,i)=u;
  end
  for j=2:mm-1
    h=HH(j-1,:);
    Ld=find(h>=0);
    v=squeeze(VV(k,j,:));
    v(end)=nan;
    v(Ld)=nan;
    Vc(k,j,:)=v;
  end
end



return