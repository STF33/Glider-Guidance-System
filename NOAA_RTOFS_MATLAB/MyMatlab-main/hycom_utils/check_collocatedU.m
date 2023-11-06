   function colloc = check_collocatedU(HH,U,V);
% colloc = check_collocatedU(HH,U,V);
%
% Input U, V should be nan on land!
%
% In HYCOM output fields U and V can be collocated or not
%  Computational grid is C:
%
%     |             |
%     |      i,j    |
%     -U(i,j) *     |
%     |       Tr    |
%     |             |
% -------- | V(i,j)-------
%
% usually, archv, archm files on original vertical HYCOM grid
% are not collocated, fields interpolated on z-levels are 
% collocated
% 
% This functino does a quick check of whether the U and V are 
% collocated
%
% u(i,j) is between h(i-1,j) and h(i,j)
% To check look at near land points
% such that HH(j-1,i) and HH(j,i-1) is land and HH(j,i) is not
% if V(j,i) is nan  than V is not collocated with H
% if U(j,i) is nan than U is not collocated with H
%
% Input: 
% HH - topo, land >0
% U, V - velocity fields, 2D, surface layer
% land values - nan's or 1e30
%
% Output:
% colloc - logical 1 - collocated, 0 - not collocated 
%
[mm,nn]=size(HH);
ndm = ndims(U);
if ndm==3
  U=squeeze(U(1,:,:));
  V=squeeze(V(1,:,:));
end


U(U>1e20)=nan;
V(V>1e20)=nan;

% Checking U:
j  = 25;
i  = 2;
h0 = -100;
while h0<0
  j=j+1;
  h0=HH(j,i);
end
j0=j;
while h0>=0
  i=i+1;
  h0=HH(j,i);
end
i0=i;

h0=HH(j0,i0);
u0=U(j0,i0);  
%keyboard
% check land value:
if ~isnan(U(j0,i0-1)),
  error('check_collocatedU: land value for U is not nan or huge');
end

if isnan(h0),
  fprintf('j0=%i, i0=%i, h0=%d\n',j0,i0,h0);
  error('Troulbe finding not-land point ...');
end

if isnan(u0)
  clc1=logical(0); % not collocated;
else
  clc1=logical(1); % collocated
end

% check V:
j  = 10;
i  = round(nn/4);
i0 = i;
h0 = -100;
while h0<0
  j=j+1;
  h0=HH(j,i);
end
while h0>=0
  j=j+1;
  h0=HH(j,i);
end
j0=j;
h0=HH(j0,i0);
v0=V(j0,i0);
 
if isnan(h0),
  fprintf('j0=%i, i0=%i, h0=%d\n',j0,i0,h0);
  error('Troulbe finding not-land point ...');
end

if isnan(v0)
  clc2=logical(0); % not collocated;
  fprintf('U,V, Not collocated\n');
else
  clc2=logical(1); % collocated
  fprintf('U,V are collocated\n');
end

if clc1~=clc2
  error('Contradiction between U and V points');
end

colloc = clc1;

return