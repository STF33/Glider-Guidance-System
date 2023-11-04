function CMP = create_colormap8(nint,c1,c2);
% CMP = create_colormap7(nint,c1,c2);
% Create colormap: very smooth transition
%
% 6 color groups x 4 within each 
% Creates colormap:  white - blue - green- yellow -dark red
%
%nint=24;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
%      0.7 0. 0.5;
sclr=[ 1 1 1;
      0.4 0.4 1;
      0.4 0.8  1;
      0.6 0.9  0.6;
      1 0.9 0.4;
      1 0.7 0.4;
      0.8 0.3 0.1];

%sclr = sclr/255;
nsb=size(sclr,1)-1;
ni=round(nint/nsb);

cmp=[];
for ik=1:nsb;
  cl1=sclr(ik,:);
  cl2=sclr(ik+1,:);
  clrM=mix_2colors(cl1,cl2,ni);
  icc=size(cmp,1);
  cmp(icc+1:icc+ni,:)=clrM;
end;

cmp = smooth_colormap(cmp,15);
cmp = smooth_colormap(cmp,15);
  
CMP.colormap=cmp;
CMP.intervals=cnt;

return
