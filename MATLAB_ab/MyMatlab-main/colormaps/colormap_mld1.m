  function CMP = colormap_mld1(nint,c1,c2);
% Create colormap: darkblue - turqouise - yellow - white

cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

% Create colormap by mixing 2 colors
% within 1 interval
%
sclr=[0, 0, 0.3;
      0.2 0.15 0.5;
      0, 0.4, 0.89;  % 
      0.05, 0.58, 0.83;
      0.11, 0.7, 0.69;
      0.53, 0.75, 0.47;
      0.89, 0.73, 0.32;
      0.98, 0.98, 0.05;
      0.7, 0.7, 0.7;
      0.9, 0.9, 0.9;
      1, 1, 1];
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
  
CMP.colormap=cmp;
CMP.intervals=cnt;
