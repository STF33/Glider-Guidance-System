% Vertically interpolate GLORYS fields into NEMO layer 
% between the two GLORYS layers
% Use basis functions: T(z) = T(z1)*phi1(z)+T(z2)*phi2(z)
% where phi = a0+a1*z, and phi1(z1)=1, phi1(z2)=0, phi2(z1)=0, phi2(z2)=1
function AA = sub_vrtintrp_glorys(A1,A2,zz1,zz2,zz0);

ZP = [1    1; ...
      zz1 zz2];
II = eye(2);
invZP = 1./(zz2-zz1)*[zz2, -1; -zz1, 1];
PHI = II*invZP;

% Check:
%t1=-1;
%t2=5;
%zm1=[1; zz1];
%b = PHI*zm1;
%t1i = b'*[t1;t2]; % should be t1
%zm2=[1; zz2];
%b = PHI*zm2;
%t2i = b'*[t1;t2];  % should be t2
%zm3 = [1; 0.5*(zz1+zz2)];
%b = PHI*zm3;
%t3i = b'*[t1;t2];  % should be mean
%keyboard

zm1 = [1;zz0];
W = PHI*zm1;    % basis fns at depth zm1, i.e. weights
AA = W(1)*A1+W(2)*A2;

return
