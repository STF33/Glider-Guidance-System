% % Apply 2D Gaussian filter
% Near NEMO/GLORYS boundaries where
% there is a jump in the 2D fields
%
% Input: nij - half-dim of the filter, 
%              i.e.# of points from the center 
%        eM1 - input field to filter
% smooth only ourside of the NEMO domain
% di points, Rflt - Gaussian "radius"
% No smoothing on the West bndry - land
%        
% Output: eMf - filtered field
function eMf=sub_fltr_Gauss(nij,eM1,jbot,jtop,ilft,irht,din,dout);

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

%keyboard

% Filter:
fprintf('NEMO boundaries, Gaussian Filtering, %ix%i pnts ...\n',ndm,ndm);
[mm,nn]=size(eM1);
ntt=mm*nn;
eMf=eM1;
cc=0;

% Bottom:
for ibnd=1:3
  if ibnd==1  % Bottom
    ii1=ilft-1;
    ii2=irht+1;
    jj1=jbot-dout;
    jj2=jbot+din;
  elseif ibnd==2; % Eastern bndry
    ii1=irht-din;
    ii2=irht+dout;
    jj1=jbot-1;
    jj2=jtop+1;
  elseif ibnd==3; % North
    ii1=ilft-1;
    ii2=irht+1;
    jj1=jtop-din;
    jj2=jtop+dout;
  end    

  for ii=ii1:ii2
    for jj=jj1:jj2
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
%keyboard
				end
		end

end

%keyboard 

return

