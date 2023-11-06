function CMP = create_colormapBGY(nint,c1,c2);
% CMP = create_colormapBGY(nint,c1,c2);
% Create colormap: very smooth transition
%
% 6 color groups x 4 within each 
% Creates colormap:  dark- blue - green- yellow 
%
%nint=24;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
%      0.7 0. 0.5;
sclr=[0 0 0.3
      0 0.2 0.3
      0 0.3 0.4
      0 0.4 0.5
      0 0.5 0.6
      0.2 0.5 0.6
      0.3 0.6 0.5
      0.5 0.7 0.4
      0.9 0.8 0.4 
      1 0.8 0.4
      1 0.8 0.5
      1 1 0.8];

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

if nint>15
cmp = smooth_colormap(cmp,15);
cmp = smooth_colormap(cmp,15);
end

CMP.colormap=cmp;
CMP.intervals=cnt;

return
