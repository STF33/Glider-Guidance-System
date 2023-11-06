 function [II,JJ]=sub_xsct_indx(is,js,ie,je);
%   [I,J]=sub_xsct_indx(i1,j1,i2,j2);
% indices - Integers
% for given end-points of a cross-section
% get indices along the xsection line
% that connectst the endpoints


%dmm1=0;
%if ie<is
%  dmm1=ie;
%  ie=is;
%  is=dmm1;
%end

clear II JJ
II(1,1)=is;
JJ(1,1)=js;
cntr=1;
icnt=1;
x0=is;
y0=js;

%keyboard

ss0=sqrt((is-ie).^2+(js-je)^2);
ss=1.5*ss0;

while ss>0.5
  ss=sqrt((x0-ie).^2+(y0-je)^2);
  if ss<0.5; break; end;
  di=(ie-x0)/ss;
  dj=(je-y0)/ss;
  Xnrm=[-dj,di];  % normal to the X-section
  ss=sqrt((x0-ie).^2+(y0-je)^2);
  di=(ie-x0)/ss;
  dj=(je-y0)/ss;
  x0=x0+di;
  y0=y0+dj;
  if round(x0)==II(cntr) & round(y0)==JJ(cntr)
    continue;
  end
  
  cntr=cntr+1;
  II(cntr,1)=round(x0);
  JJ(cntr,1)=round(y0);
  if (cntr>1e6); 
    fprintf('Couldnt find indices? Too many points: %s\n',cntr);
    keyboard;
%    break; 
  end;
  if ss>2*ss0
    error('Points are moving away from the segment');
  end
  
end

%keyboard

%if dmm1>0
%  II=flipud(II);
%  JJ=flipud(JJ);
%end

% Make sure that di<=1 and dj<=1
% note: works only in index space
f_didj=1;
ni=length(II);
cc=0;
if f_didj
for k=1:ni-1;
  i1=II(k);
  i2=II(k+1);
  j1=JJ(k);
  j2=JJ(k+1);
  
  nl=sqrt((i2-i1).^2+(j2-j1).^2);
  cc=cc+1;
  IIo(cc)=II(k);
  JJo(cc)=JJ(k);
  if nl<=1, continue; end
  di=i2-i1;
  dj=j2-j1;
  Np=abs(di)+abs(dj)-1;
  iold=i1;
  jold=j1;
%  if abs(di)>1 | abs(dj)>1; pause; end;
  
  for mp=1:abs(di)
    cc=cc+1;
    IIo(cc)=IIo(cc-1)+sign(di);
    JJo(cc)=JJo(cc-1);
  end

  for mp=1:abs(dj)
    cc=cc+1;
    IIo(cc)=IIo(cc-1);
    JJo(cc)=JJo(cc-1)+sign(dj);
  end
% Should arrive at pnt k+1:
  nl=sqrt((IIo(cc)-i2).^2+(JJo(cc)-j2).^2);
  if nl>0
    error('Missed the final point');
  end
  cc=cc-1;
end  
cc=cc+1;
IIo(cc)=II(end);
JJo(cc)=JJ(end);
II=IIo;
JJ=JJo;
end
 
 
 return
