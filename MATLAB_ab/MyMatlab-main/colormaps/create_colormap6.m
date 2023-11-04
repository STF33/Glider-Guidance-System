  function CMP = create_colormap6(nint,c1,c2);
% CMP = create_colormap5(nint,c1,c2);
% Creates colormap: blue - green- yellow -red
% 
%  example: 
% nint=200;
%     c1=0;
%     c2=20;
%  CMP=create_colormap5(nint,c1,c2);
%  cmp=CMP.colormap;
%  cnt=CMP.intervals;


%nint=100;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

% Create colormap by mixing 2 colors
% within 1 interval
%
sclr=[0 0.4 0.7;
      0 1   1;
      0 1  0;
      1 1 0;
      0.5 0.5 0;
      1 0 0;
      0.4 0 0];

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
  
CMP.colormap=cmp;
CMP.intervals=cnt;
