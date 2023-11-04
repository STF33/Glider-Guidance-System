  function CMP = create_PBYR(nint,c1,c2);
% CMP = create_PBYR(nint,c1,c2);
% Create colormap: purple-blue-yellow-red
%
%  CMP=create_colormap1(nint,c1,c2);
%  cmp=CMP.colormap;
%  cnt=CMP.intervals;


%nint=100;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

% Create colormap by mixing 2 colors
% within 1 interval
sclr=[0.7 0 0.7;
      1   0  1;
      1   0.8 1;
      0.8 0.8 1;
      0   0.5 1;
      0   0.3 0.8;
      1, 1, 0.8;
      1, 0.8, 0.4;
      1, 0.6, 0;
      1, 0.8, 0.8; 
      1, 0.4  0.4;
      0.8, 0, 0];

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

cmp=smooth_colormap(cmp,15);
cmp=smooth_colormap(cmp,15);
  
CMP.colormap=cmp;
CMP.intervals=cnt;
