function CLC=sub_collocate2uv(U,T,S,dH,HH,i0,j0,i1,j1);
% collocate p-point variables with U or V pnt
% only 2D and 3D allowed
dm=size(U);
ndm=max(size(dm));
if ndm>3, error('sub_collocate2uv.m: Dimension > 3'); end
%keyboard
if ndm==3
  un=U(:,j0,i0);
  t1=T(:,j0,i0);
  t2=T(:,j1,i1);
  s1=S(:,j0,i0);
  s2=S(:,j1,i1);
  dh1=dH(:,j0,i0);
  dh2=dH(:,j1,i1);
else
  un=U(j0,i0);
  t1=T(j0,i0);
  t2=T(j1,i1);
  s1=S(j0,i0);
  s2=S(j1,i1);
  dh1=dH(j0,i0);
  dh2=dH(j1,i1);
end
h1=HH(j0,i0);
h2=HH(j1,i1);

% Near land points:
if h2>=0
  t2=t1;
  s2=s1;
end


CLC.Un=un;
CLC.Tn=0.5*(t1+t2);
CLC.Sn=0.5*(s1+s2);
dHn=0.5*(dh1+dh2);
sH=nansum(dHn);
Hn=abs(0.5*(h1+h2));
CLC.dHn=dHn;
CLC.Hn=Hn;

if abs(sH)>0.01 & abs(Hn-sH)/sH>0.01
fprintf('** ERR sub_collocate2uv:  depths i=%i, j=%i, H=%6.2f dh=%6.2f\n',...
	    i0,j0,Hn,sH);
%    keyboard
end



return