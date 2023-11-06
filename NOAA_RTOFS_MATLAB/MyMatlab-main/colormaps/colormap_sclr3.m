  function CMP = colormap_sclr3(nint,c1,c2);
% Create colormap: jet with light pink at the end

cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

% Create colormap by mixing 2 colors
% within 1 interval
%
sclr =[0    0      0.5
       0    0.5    1
       0    1      1
       0.5  1      0.5
       1    1      0
       1    0.5    0
       1    0      0
       1.0  2      2];

nsb=size(sclr,1)-1;
ni=ceil(nint/nsb);
%keyboard

cmp=[];
for ik=1:nsb;
  cl1=sclr(ik,:);
  cl2=sclr(ik+1,:);
  clrM=mix_2colors(cl1,cl2,ni);
  icc=size(cmp,1);
  cmp(icc+1:icc+ni,:)=clrM;
end;
  
cmp = smooth_colormap(cmp,ni);

CMP.colormap=cmp;
CMP.intervals=cnt;

return