% Find color for given value
function CC = sub_match_color(ssh,CMP,clm1,clm2);

ncmp = size(CMP,1);
% intervals:
dx = (clm2-clm1)/ncmp;
xint = [clm1:dx:clm2];


ns = length(ssh);
for ik=1:ns
  s0 = ssh(ik);
  if isnan(s0); 
    CC(ik,:)=[NaN,NaN,NaN];
    continue;
  elseif s0<=clm1; 
    CC(ik,:)=CMP(1,:);
    continue;
  elseif s0>=clm2;
    CC(ik,:)=CMP(end,:);
    continue;
  end

  jx = max(find(xint<=s0));
  CC(ik,:)=CMP(jx,:);
end

return
