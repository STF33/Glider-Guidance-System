function CMP = colormap_discrete01(c1,c2);
% Create discrete colormap
% 6 color groups x 4 within each 
% Creates colormap: blue - green- yellow -red
%
nint=24;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
%      0.7 0. 0.5;
sclr=[ 1 1 1;
      0. 0. 1;
      0 0.5   0.6;
      0 0.7  0;
      0.6 0.6 0;
      1 0.5 0;
      0.4 0 0];

%sclr = sclr/255;
nsb=size(sclr,1)-1;
ni=round(nint/nsb);

cmp=[];
for ik=1:nsb;
%  cl1=sclr(ik,:);
  cl1=[1 1 1];
  cl2=sclr(ik+1,:);
  clrM=mix_2colors(cl1,cl2,ni);
  clrM=smooth_colormap(clrM,3);
  clrM=smooth_colormap(clrM,3);
  icc=size(cmp,1);
  cmp(icc+1:icc+ni,:)=clrM;
end;
  
CMP.colormap=cmp;
CMP.intervals=cnt;

return