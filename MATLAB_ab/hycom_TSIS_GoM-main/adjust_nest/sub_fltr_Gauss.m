function eMf=sub_fltr_Gauss(nij,eM1);
% Apply 2D Gaussian filter
% Near boundaries - filter is 
% choped (i.e. non symmetric)
% however the overall weight is kept = 1
%
% Input: nij - hlaf-dim of the filter, 
%              i.e.# of points from the center 
%        eM1 - input field to filter
%        
% Output: eMf - filtered field
ndm=nij*2+1;
sgm=nij;
sgm2=sgm^2;
SGM2=zeros([ndm,1])+sgm^2;
SGM2=diag(SGM2,0);
ii=0;
clear GS0
for i=-nij:nij
  ii=ii+1;
  jj=0;
  for j=-nij:nij
    jj=jj+1;
    x0=sqrt(i^2+j^2);
    expn=-(x0^2)/(2*sgm2);
    GS0(jj,ii)=exp(expn);
%    fprintf('jj=%i, ii=%i, x0=%6.4f, GS=%6.4f\n',jj,ii,x0,exp(expn));
  end
end

GS=GS0./sum(sum(GS0));
[mf,nf]=size(GS);
%pcolor(GS); shading flat
%colorbar

% Filter:
fprintf('Gaussian Filtering, %ix%i pnts ...\n',ndm,ndm);
[mm,nn]=size(eM1);
ntt=mm*nn;
eMf=eM1*nan;
cc=0;
for ii=1:nn
  for jj=1:mm
    cc=cc+1;
    if mod(cc,150000)==0,
      fprintf('::: Done %5.1f%% ...\n',cc/ntt*100);
    end
    
    if isnan(eM1(jj,ii)), continue; end;
    i1=ii-nij;
    i2=ii+nij;
    j1=jj-nij;
    j2=jj+nij;
    if1=1;
    if2=nf;
    jf1=1;
    jf2=mf;
    if i1<1
      di=-i1+1;
      if1=1+di;
    end
    if i2>nn
      di=i2-nn;
      if2=nf-di;
    end
    if j1<1
      dj=-j1+1;
      jf1=1+dj;
    end
    if j2>mm
      dj=j2-mm;
      jf2=mf-dj;
    end
    A=GS(jf1:jf2,if1:if2);
    A=A./sum(sum(A));
    
    i1=max([1,i1]);
    i2=min([i2,nn]);
    j1=max([1,j1]);
    j2=min([j2,mm]);
    
    B=eM1(j1:j2,i1:i2);
    dmm=nansum(nansum(B.*A));
    eMf(jj,ii)=dmm;  
    
  end
end
I=find(isnan(eM1));
eMf(I)=nan;
fprintf('Done filtering\n');


return