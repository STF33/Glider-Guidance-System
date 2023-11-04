function Pn = spline_cub(Xi,Fi,Xnew,nbcs);
% Given mesh point Xi and values Fi
% fit a piecewise polynomial of degree dn
% and interpolate original data
% onto a higher-resolution mesh Xnew
% Output Pns - function at Xnew points
%
% Dmitry Dukhovskoy, FSU, 2018-2020
%
% Cubic spline
% pi(x) = s2(i-1)/6h(i)*(xi-x)^3+s2(i)/6h(i)+...
%         gamm(i-1)*(x-x(i-1))+gamt(i-1)
% gamt = ...
% gamm = ...
% Need to find coefficients, s2(0,...,n)
% For this, solve T*s2=d
%
% + two BCs

% Boundary conditions:
if isempty(nbcs)
  nbcs = 'natural'; % 
end
%nbcs = 'herm1'; % s0'=f'(x0), sn'=f'(xn),Hermite bcs1 
%nbcs = 'herm2'; % s0"=f"(x0), sn"=f"(xn), Hermite bcs2


np = length(Xi); % total # of points = n+1, Xi(1) = X0
dn = 3; % set the degree of the (cubic) spline

% Create points where to evaluate the polynomial:
%dx = 0.01;
%xx = [x0:dx:xN];
xx=Xnew;
% Define H step
for ik=1:np
  if ik==1
    H(ik) = Xi(ik+1)-Xi(ik);
  else
    H(ik) = Xi(ik)-Xi(ik-1);
  end
end  


%[s2_0,s2_n] = sub_bcs.m % get s"(0) and s"(n) for different BCs
switch(nbcs)
 case('natural');
  dmus0=0; % OBc
  dmusN=0;
 case('herm1'); % s'(0)=f'(0), s'(n)=f'(n)
  ik=1;
  mu  = H(ik)/(H(ik)+H(ik+1));
  a = (Xi(1)-Xi(2)).^2./(2*H(1));
  b = (Fi(2)-Fi(1))/H(2);
  c = H(2)/6;
  cff0 = 1/(1+c/(2*(c-a))); % convert aprx S" -> s"
  dmus0 = -mu*(f1_0-b)/(c-a);

  ik=np-1;
  lmb = H(ik+1)/(H(ik)+H(ik+1));
  a = (Xi(np)-Xi(np-1)).^2./(2*H(np));
  b = (Fi(np)-Fi(np-1))/H(np);
  c = H(np)/6;
  cffn = 1/(1-lmb*c/(2*(a-c)));
  dmusN = lmb*(b-f1_n)/(a-c);
  
  
 case('herm2');
  ik=1;
  mu  = H(ik)/(H(ik)+H(ik+1));
  dmus0 = -mu*f2_0;
  ik=np-1;
  lmb = H(ik+1)/(H(ik)+H(ik+1));
  dmusN = -lmb*f2_n;
end


% Natural BCs
% T matrix is (n-1)x(n-1) 
% s2(0)=s2(n) = 0

T = zeros(np-2,np-2); % gives (n-1)x(n-1) matrix, where 0<i<n x0=x(1)=0, etc.
S = zeros(np-2,1);
D = zeros(np-2,1);
dmT = size(T,1);
cc = 0;
for ik=2:np-1
  cc = cc+1;
  
  lmb = H(ik+1)/(H(ik)+H(ik+1));
  mu  = H(ik)/(H(ik)+H(ik+1));
  dd  = 6/(H(ik)+H(ik+1))*((Fi(ik+1)-Fi(ik))...
			   /H(ik+1)-(Fi(ik)-Fi(ik-1))/H(ik));

  if ik==2,
    dd = dd+dmus0; % effect of the BCs at i=0
  end
  if ik==np-2
    dd = dd+dmusN; % BC at i=n
  end
  
  
  imu = cc-1;
  i2 = cc;
  ilmb = cc+1;
  if imu>0
    T(cc,imu) = mu;
  end;
  if i2<=dmT, T(cc,i2)=2; end;
  if ilmb<=dmT, T(cc,ilmb)=lmb; end;
  
  D(cc,1) = dd;
end

%keyboard

% Solve for S
% Do LU
%[L,U] = lu(T);

% Find 
S2 = inv(T)*D;

if strcmp(nbcs,'herm1'), % adjust s1 and s(n-1)
  S2(1) = S2(1)*cff0;
  S2(end) = S2(end)*cffn;
end


%
% Add missing rows to have correct dims
% for natural BC
S2=[0;S2;0];
clear Pn
Pn = [];
xns = []; % bookeeping, xns should = xx after all done
%Pn(1) = Fi(1); % sontraint at the OB
%Pn(np) = Fi(np);
% pi(x) = A*(x-xi^3)+B*(x-xi-1^3)+gamm*(x-xi-1)+gamt;
for ik=2:np
  A = S2(ik-1)./(6*H(ik));
  B = S2(ik)./(6*H(ik));
  gamt = Fi(ik-1)-S2(ik-1)*H(ik)^2/6;
  gamm = (Fi(ik)-Fi(ik-1))/H(ik)-H(ik)/6*(S2(ik)-S2(ik-1));
  
% Find points inside the segments
  x0s = Xi(ik-1);
  xNs = Xi(ik);
  ix1 = min(find(xx>=x0s & xx<xNs));
  ix2 = max(find(xx<=xNs & xx>x0s));
  if isempty(ix1) | isempty(ix2)
    fprintf('Segment %i contains no data points\n',isgm);
    continue
  end
  xxs = xx(ix1:ix2); % subsample xx points for this segment
  xxs = xxs(:);

  Pi = A*(Xi(ik)-xxs).^3+B*(xxs-Xi(ik-1)).^3+...
       gamm*(xxs-Xi(ik-1))+gamt;
  Pi = Pi(:);
  
  if ~isempty(xns) & xxs(1)==xns(end);
    xxs=xxs(2:end);
    Pi=Pi(2:end);
  end
  
  Pn = [Pn;Pi];
    
  xns = [xns;xxs];
%  keyboard
end


f_chck=0;
if f_chck
figure(1); clf;
plot(Xi,Fi,'o-');
hold on;
plot(xns,Pn,'r');
stt = sprintf('Cubix Spline (red)');
title(stt);
end


return
