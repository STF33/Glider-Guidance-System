function Fdf = sub_newton_bern(Xi,Fi);
% Newton polynomila for a set of points:
% [x0,...,xn], on [a,b], x0=a, xn=b
% Use Smoktunowicz alg. to compute 
% divided difference algorithm
% this is O(n) computations and space
% Input Xi and Fi  - 1D arrays
% 
% Dmitry Dukhovskoy, FSU, 2018
%
np = length(Fi);
dn = np-1; % degree of the polynomial
% Create array of f[x0,x1,...,xn] 
% divided differences
% This method gives O(n^2) operations
% To create an array of f[x0,...,xn]
% to update from f[x0,...,xn-1] -> f[x0,...,xn] - O(n) operations
Fdf = zeros(np,1);
Fdf(1) = Fi(1); % f0
cc = 0;
for kk=2:np % loop for differences, kk=2 ->f[x0,x1],...
  ff = 0;
  for ii=1:kk % loop for omg_kk(x(ii))
% calculating omg_kk(x(ii))    
    omg = 1;
    for ip=1:kk % loop kk=3,ii=1, ip=1:3 =>
                % omg_3(x0)=(x0-x1)*(x0-x2), (x0-x0)- skipped
      if ii==ip, continue; end; % skip ii=kk in omg calc
      dmm=(Xi(ii)-Xi(ip));
      omg=omg*dmm;
      cc = cc+1;  % counter
    end
    ff =ff+Fi(ii)/omg; % f[x0,...,x(kk-1)]
%    keyboard
  end
  Fdf(kk)=ff;
end
%Fdf2 = Fdf;

%keyboard

return
