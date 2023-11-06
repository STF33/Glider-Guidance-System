function [Uc,Vc] = sub_collocateUV2p_1lr(U,V,dH,HH);
% collocate U and V components
% into p-points
% 1 layer only
%
fprintf('sub_collocateUV2p_1lr: Collocating UV to p-points\n');
dm=size(U);
ndm=max(size(dm));
if ndm>2, error('Dimension > 2'); end

[mm,nn]=size(U);
for jj=2:mm-1
  h1=HH(jj,:);
  h2=HH(jj+1,:);

  I=find(h1>=0 & h2<0);
  h1(I)=h2(I);

  I=find(h2>=0 & h1<0);
  h2(I)=h1(I);

  u1=U(jj,:);
  u2=U(jj+1,:);

  dh1 = 0.5*(dH(jj-1,:)+dH(jj,:));
  dh2 = 0.5*(dH(jj,:)+dH(jj+1,:));
  
  Uc(jj,:)=(dh1.*u1+dh2.*u2)./(dh1+dh2);
end
Uc(1,:)=U(1,:);
Uc(mm,:)=U(mm,:);

for ii=2:nn-1
  h1=HH(:,ii);
  h2=HH(:,ii+1);

  I=find(h1>=0 & h2<0);
  h1(I)=h2(I);

  I=find(h2>=0 & h1<0);
  h2(I)=h1(I);

  v1=V(:,ii);
  v2=V(:,ii+1);

  dh1 = 0.5*(dH(:,ii-1)+dH(:,ii));
  dh2 = 0.5*(dH(:,ii)+dH(:,ii+1));

  Vc(:,ii)=(dh1.*v1+dh2.*v2)./(dh1+dh2);
end
Vc(:,1)=V(:,1);
Vc(:,nn)=V(:,nn);



return
