    function [UA,IP] = check_spuriousU(dmm);
%
% In HYCOM mean fields, spurious velocities can exist
% need to check
% The erroneous points are identified and filled with values from neighbouring points
%  dmm - U or V field with NaNs for land mask
%   


[m,n]=size(dmm);
FX = diff(dmm,1,2);
FX = [zeros(m,1),FX];
FY = diff(dmm,1,1);
FY = [zeros(1,n);FY];

S = sqrt(FX.*FX+FY.*FY);
[m,n]=size(dmm);

% Delete the most exterme values
IK = find(S>0.35);

umm = dmm;
umm(IK)=nan;
umm(umm==0) = nan;  % layers with 0 thicknesses

IK = find(S>0.35);
[J,I] = ind2sub(size(S),IK);

UA = dmm;
d=3;
ck=0;
IP=[];
for ii=1:length(IK)
  i0=I(ii);
  j0=J(ii);
  gr=S(j0,i0);
  u=dmm(j0,i0);
% Find spots around with true velocities and check if 
% suspicious point is spurious value or realistic
  im1 = max([1, i0-d]);
  ip1 = min([n, i0+d]);
  jm1 = max([1, j0-d]);
  jp1 = min([m, j0+d]); 
  T = abs(umm(jm1:jp1,im1:ip1));
  [mt,nt]=size(T);
%  P = dmm(jm1:jp1,im1:ip1);

  sb = S(jm1:jp1,im1:ip1);

  mnU = nanmean(reshape(T,mt*nt,1));
  stdU = nanstd(reshape(T,mt*nt,1));

  if abs(u)>(mnU+5*stdU);
    ck=ck+1;
    IP(ck)=IK(ii);
    un = sign(u)*mnU;
%    fprintf('Adj. spurious U: %i i0=%i j0=%i, Old u = %9.3f ===> unew = %9.3f\n',ck, i0,j0,u,un);
    UA(j0,i0)=un;
  end

end


