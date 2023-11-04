  function rgm = gamma_fn(alf0,lambda);
% Fast method for 
% generating gamma Random variable
% Marsaglia and Tsang, A simple method for
% generating Gamma variables
% Transactions on Mathematical Software, 26(3), 363-372
% adopted for alfa>1 and alfa<1
% 
% See also: 
% http://www.hongliangjie.com/2012/12/19/how-to-generate-gamma-random-variables/
% Note lambda here is 1/betta used in Matlab:
% betta=sgm/sqrt(alf); % for 1/(alf*[betta^alf]*Gamma)* ... 
% lambda=1/betta;      % for [lmbd^alf]/(alf*Gamma)* ...
%
% Example:
%dmean=20; % micrometer - mean size
%cv=0.6;   % coefficient of variation - related to alfa and mean
%alf = (1/cv)^2;
%sgm=dmean*cv; % st. dev. 
%betta=sgm/sqrt(alf); % for 1/(alf*[betta^alf]*Gamma)* ... 
%lambda=1/betta;      % for [lmbd^alf]/(alf*Gamma)* ...
%for ik=1:50000
%  rgm=gamma_fn(alf,lambda);
%  RR(ik)=rgm;
%end


if alf0<=1
  alf=alf0+1;
else
  alf=alf0;
end

rgm=-100;
cc=0;
while rgm<0
  cc=cc+1;
  if cc>10,
    error('gamma_fn: failed to generate Gamma rv for alf=%8.4f',alf);
  end
  
%  alf=11;
  d=alf-1/3;
  c=1/sqrt(9*d);
  xn=-100/c;
  while xn<=-1/c
    xn=randn;  % normal
  end
  v=(1+c*xn)^3;
  uu=rand; % uniform

%  dmm=1-0.0331*xn^4;
%  if uu<dmm
%    rgm=d*v/lambda;
%    return;
%  end

  lmm=0.5*xn^2+d*(1-v+log(v));
  if log(uu)<lmm
    rgm=d*v/lambda;
    
    if alf0<=1
      rgm=rgm*rand^(1/alf0);
    end
    
    return;
  end

end

return





