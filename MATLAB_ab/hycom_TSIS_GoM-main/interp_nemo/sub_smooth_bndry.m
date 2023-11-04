% NEMO - GLORYS boundary has a jump
% fit a poynomial to smooth the 
% discontinuity
% 
% nij - how many points  to smooth with an interpolant
% din - # of points inside the domain, 0 - 1 pnt outside domain
% jbot, jtop, ilft, irht - boundary indices
%
function eMf=sub_smooth_bndry(nij,din,eM1,jbot,jtop,ilft,irht);

fprintf('Smoothing NEMO bndry jumps ...\n');
eMf = eM1;

if nij<din, error('# of insde points > total # of points'); end;
nhlf = round(nij/2);

di = 5; % pnts outside interpolation
% Bottom:
% fix pnt din pnts in the domain, =0 - on the bndry
ii1=ilft;
ii2=irht;
for ii=ilft:irht
  jj=jbot;
  if isnan(eM1(jj,ii)), continue; end;
% Find start-end points for interpolation
  jIE = jbot+din;
  jIS = jIE-nij+1;
% if around bndry land - also skip
  I=find(isnan(eM1(jIS:jIE,ii)));
  if ~isempty(I); continue; end; 
  
  j1=jIS-di;
  j2=jIE+di;

% whole section:
  dmm = eM1(j1:j2,ii);
  xx1 = [j1:j2]';
% -----------

  xx = [j1:jIS-1, jIE+1:j2]';
%  amm = [eM1(j1:jin-1,ii); eM1(jin+nij:j2,ii)];
  amm = eM1(xx,ii);
  xxi = [jIS:jIE]';
  I=find(isnan(amm));
  if ~isempty(I); continue; end;

  ammi = interp1(xx,amm,xxi,'spline');
  eMf(xxi,ii)=ammi;
end;

%
% Eastern bndry
% 1 pnt outside
jj1=jbot;
jj2=jtop;

for jj=jj1:jj2
  ii=irht;
  if isnan(eM1(jj,ii)), continue; end;
  jIS = irht-din;
  jIE = irht+nij;
% if around bndry land - also skip
  I=find(isnan(eM1(jIS:jIE,ii)));
  if ~isempty(I); continue; end;

  i1=jIS-di;
  i2=jIE+di;
%keyboard
  xx = [i1:jIS-1,jIE+1:i2]';
  amm = eM1(jj,xx)';
  xxi = [jIS:jIE]';
  I=find(isnan(amm));
  if ~isempty(I); continue; end;

  ammi = interp1(xx,amm,xxi,'spline');
  eMf(jj,xxi)=ammi';

%  dmm = eM1(jj,i1:i2);
%  xx1 = [i1:i2]';

end;

%
% North bndry
% din = 0 - starts 1 pnt outside of bndry
for ii=ilft:irht
  jj=jtop;
  if isnan(eM1(jj,ii)), continue; end;
  jIS = jtop-din;
  jIE = jIS+nij;
% if around bndry land - also skip
  I=find(isnan(eM1(jIS:jIE,ii)));
  if ~isempty(I); continue; end;

  j1=jIS-di;
  j2=jIE+di;

  xx = [j1:jIS-1, jIE+1:j2]';
%  amm = [eM1(j1:jin-1,ii); eM1(jin+nij:j2,ii)];
  amm = eM1(xx,ii);
  xxi = [jIS:jIE]';
  I=find(isnan(amm));
  if ~isempty(I); continue; end;

  ammi = interp1(xx,amm,xxi,'spline');
  eMf(xxi,ii)=ammi;

%  dmm = eM1(j1:j2,ii);
%  xx1 = [j1:j2]';

end;

  

return
